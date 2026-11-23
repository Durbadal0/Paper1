// sbart_engine.cpp -- C++ version of Steps 1-3 of the simulation worker.
//
// sbart_step1() is the Step-1 loop of sim_worker.R (CAR model for M) and
// sbart_step2() is its Step-2/Step-3 loop (SBART survival model, importance weights,
// AMSE and AES running sums). Both draw their random numbers from R's generator,
// in the same order as the R loops, and compute every quantity in the same order
// of operations, so that sim_worker_cpp.R reproduces sim_worker.R.
//
// The large prediction blocks (AMSE grid, AES grid, Step-3 quadrature) are computed
// on several threads. Each row is computed exactly as Forest::do_predict computes
// it, so the number of threads does not change any result.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Random.h>
#include <thread>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>

#include "softbart/soft_bart.h"
#include "softbart/soft_bart_impl.h"

// ---------------------------------------------------------------------------
// Thread-safe prediction. Forest::do_predict stores the leaf weights in the nodes,
// so it cannot run on several threads; tree_acc() keeps them on the stack. The
// arithmetic is the same: left weight = w * cw, right weight = (1 - w) * cw, leaf
// contributions added left to right, trees added in order, both sums starting at 0.
static void tree_acc(const Node* n, const double* x, double cw, double& acc) {
  if (n->is_leaf) {
    acc = acc + cw * n->mu;
    return;
  }
  double weight = activation(x[n->var], n->val, n->tau);
  tree_acc(n->left, x, weight * cw, acc);
  tree_acc(n->right, x, (1 - weight) * cw, acc);
}

static void predict_range(const std::vector<Node*>& trees, const double* X, int nrow,
                          int ncol, int r0, int r1, double* out) {
  std::vector<double> xr(ncol);
  int num_tree = trees.size();
  for (int i = r0; i < r1; i++) {
    for (int j = 0; j < ncol; j++) xr[j] = X[i + (size_t) nrow * j];
    double tot = 0.0;
    for (int t = 0; t < num_tree; t++) {
      double acc = 0.0;
      tree_acc(trees[t], xr.data(), trees[t]->current_weight, acc);
      tot = tot + acc;
    }
    out[i] = tot;
  }
}

static void predict_mt(const Forest& forest, const double* X, int nrow, int ncol,
                       double* out, int nthreads) {
  const std::vector<Node*>& trees = forest.trees_ref();
  if (nthreads <= 1 || nrow < 2000) {
    predict_range(trees, X, nrow, ncol, 0, nrow, out);
    return;
  }
  std::vector<std::thread> pool;
  int chunk = (nrow + nthreads - 1) / nthreads;
  for (int k = 0; k < nthreads; k++) {
    int r0 = k * chunk, r1 = std::min(nrow, r0 + chunk);
    if (r0 >= r1) break;
    pool.emplace_back(predict_range, std::cref(trees), X, nrow, ncol, r0, r1, out);
  }
  for (auto& th : pool) th.join();
}

// ---------------------------------------------------------------------------
// From here on, a*b + c must be rounded as two operations, as R does it.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// truncnorm 1.0-9, rtruncnorm() for one-sided truncation at 0 (sd = 1); the same
// draws as truncnorm::rtruncnorm (checked by verify/check_truncnorm.R)
static const double TN_T4 = 0.45;
static inline double tn_ers_a_inf(double a) {
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    x = R::rexp(ainv) + a;
    rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0, 1) > rho);
  return x;
}
static inline double tn_nrs_a_inf(double a) {
  double x = -DBL_MAX;
  while (x < a) x = R::rnorm(0, 1);
  return x;
}
static inline double tn_left(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  if (alpha < TN_T4) return mean + sd * tn_nrs_a_inf(alpha);
  return mean + sd * tn_ers_a_inf(alpha);
}
static inline double tn_right(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  return mean - sd * tn_left(-beta, 0.0, 1.0);
}

// for verify/check_truncnorm.R: draws from the sampler above
// [[Rcpp::export]]
Rcpp::NumericVector engine_rtruncnorm(Rcpp::NumericVector mean, int lower) {
  Rcpp::NumericVector out(mean.size());
  for (int i = 0; i < mean.size(); i++)
    out[i] = lower ? tn_left(0.0, mean[i], 1.0) : tn_right(0.0, mean[i], 1.0);
  return out;
}

// R's sum(), cumsum() and rowMeans() accumulate in long double
typedef long double LDOUBLE;

// car_ld() of sim_worker.R. The log-determinant term is computed by R itself
// (logdet_fun, the same eigen() call as in car_ld) and cached by the value of rho;
// the quadratic form follows the order of R's two dgemv calls for t(x) %*% Q %*% x.
struct CarLD {
  int N;
  const double* A;
  const double* D;
  Rcpp::Function logdet_fun;
  std::map<double, double> cache;
  std::vector<double> Q, z;
  CarLD(int N_, const double* A_, const double* D_, Rcpp::Function f)
    : N(N_), A(A_), D(D_), logdet_fun(f), Q(N_ * N_), z(N_) {}
  double logdet(double rho) {
    auto it = cache.find(rho);
    if (it != cache.end()) return it->second;
    double v = Rcpp::as<double>(logdet_fun(rho));
    cache[rho] = v;
    return v;
  }
  double quadform(const double* x, double rho) {
    for (int k = 0; k < N * N; k++) {
      double ra = rho * A[k];
      Q[k] = D[k] - ra;
    }
    for (int j = 0; j < N; j++) {          // z = t(Q) x   (dgemv "T")
      double temp = 0.0;
      for (int i = 0; i < N; i++) { double pr = Q[i + N * j] * x[i]; temp = temp + pr; }
      z[j] = temp;
    }
    double y = 0.0;                         // y = z x      (dgemv "N")
    for (int j = 0; j < N; j++) { double pr = x[j] * z[j]; y = y + pr; }
    return y;
  }
  double operator()(const double* x, double rho, double s2) {
    double a = logdet(rho) - N * log(s2);
    double t1 = 0.5 * a;
    double t2 = (0.5 * quadform(x, rho)) / s2;
    return t1 - t2;
  }
};

// cumtrap() of sim_worker.R for a matrix P (nr x nc, column-major), written to
// out (nr x (nc - 1)): cumulative trapezoid sums times h
static void cumtrap(const double* P, int nr, int nc, double h, double* out) {
  for (int r = 0; r < nr; r++) {
    LDOUBLE acc = 0.0;
    for (int j = 0; j < nc - 1; j++) {
      double m = (P[r + (size_t) nr * (j + 1)] + P[r + (size_t) nr * j]) / 2;
      acc += m;
      double cs = (double) acc;
      out[r + (size_t) nr * j] = cs * h;
    }
  }
}

static inline double clip01(double v) { if (v < 0) v = 0; if (v > 1) v = 1; return v; }

// ---------------------------------------------------------------------------
// Step 1: the CAR model for M, given the survey counts m_0 out of n_0i.
// [[Rcpp::export]]
Rcpp::List sbart_step1(Rcpp::List inp, Rcpp::Function logdet_fun, Rcpp::Function qf_fun) {
  Rcpp::NumericMatrix A = inp["A"], D = inp["D"];
  Rcpp::IntegerVector m0 = inp["m_0"];
  Rcpp::NumericVector n0 = inp["n_0i"];
  Rcpp::NumericVector logitM = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(inp["logitM"]));
  Rcpp::NumericVector rb = inp["rho_bounds"];
  double s0 = inp["s0"], rho0 = inp["rho0"];
  int n1 = inp["n1"];
  int N = A.nrow();
  CarLD car_ld(N, A.begin(), D.begin(), logdet_fun);

  Rcpp::NumericMatrix Msave(n1, N);
  std::vector<double> prop(N);
  for (int it = 0; it < n1; it++) {
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < N; k++) prop[k] = logitM[k];
      prop[i] = R::rnorm(logitM[i], 0.3);
      double lp = car_ld(prop.data(), rho0, s0) +
                  R::dbinom(m0[i], n0[i], R::plogis(prop[i], 0.0, 1.0, 1, 0), 1);
      double lc = car_ld(logitM.begin(), rho0, s0) +
                  R::dbinom(m0[i], n0[i], R::plogis(logitM[i], 0.0, 1.0, 1, 0), 1);
      double la = lp - lc;
      if (R_FINITE(la) && log(R::runif(0, 1)) < la)
        for (int k = 0; k < N; k++) logitM[k] = prop[k];
    }
    // qf sets the sigma^2 draw, so it is taken from R itself (R's BLAS may round the
    // sums differently from the loop in CarLD::quadform, which only enters MH ratios)
    double qf = Rcpp::as<double>(qf_fun(logitM, rho0));
    double rate = 1 + qf / 2;
    s0 = 1 / R::rgamma(1 + N / 2.0, 1 / rate);
    double rp = R::runif(std::max(rb[0], rho0 - 0.1), std::min(rb[1], rho0 + 0.1));
    double u = R::runif(0, 1);
    if (log(u) < car_ld(logitM.begin(), rp, s0) - car_ld(logitM.begin(), rho0, s0)) rho0 = rp;
    for (int k = 0; k < N; k++) Msave(it, k) = R::plogis(logitM[k], 0.0, 1.0, 1, 0);
    if ((it + 1) % 1000 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("Msave") = Msave, Rcpp::Named("logitM") = logitM,
                            Rcpp::Named("s0") = s0, Rcpp::Named("rho0") = rho0);
}

// ---------------------------------------------------------------------------
// Steps 2-3: the SBART survival model with Poisson-thinning augmentation, the
// importance weights, and the running sums for AMSE and AES.
// [[Rcpp::export]]
Rcpp::List sbart_step2(Rcpp::List inp, Rcpp::List hypers, Rcpp::List opts,
                       Rcpp::Function logdet_fun, Rcpp::Function qf_fun,
                       Rcpp::Function log_iter) {
  // ---- data
  Rcpp::NumericMatrix Xr = inp["X"], Xdr = inp["Xd"], A = inp["A"], D = inp["D"];
  Rcpp::NumericVector time = inp["time"], M_hat = inp["M_hat"], rb = inp["rho_bounds"];
  Rcpp::NumericVector sy = inp["sy"];
  Rcpp::IntegerVector delta = inp["delta"], cty = inp["cty"], d_cnt = inp["d_cnt"];
  double max_time = inp["max_time"];
  int nt = Xr.nrow(), p = Xr.ncol(), N = A.nrow(), pd = Xdr.ncol();
  // ---- settings
  int num_iter = inp["num_iter"], burn_in = inp["burn_in"], THIN = inp["THIN"];
  int use_step3 = inp["use_step3"], n_keep_max = inp["n_keep_max"];
  int nthreads = inp["nthreads"];
  // ---- prediction blocks
  Rcpp::List amse_blocks = inp["amse_blocks"];
  Rcpp::NumericMatrix aes_block = inp["aes_block"], aesM_base = inp["aesM_base"];
  int nXP = inp["nXP"], nT0 = inp["nT0"], nT = nT0 - 1, nPT0 = inp["nPT0"], nPT = nPT0 - 1;
  double dt = inp["dt"], dtp = inp["dtp"];
  Rcpp::NumericMatrix M_post = inp["M_post"];
  int NQ = inp["NQ"];
  Rcpp::NumericMatrix Xq_base = inp["Xq_base"], Xy_base = inp["Xy_base"];
  Rcpp::IntegerVector cty_q = inp["cty_q"];
  // ---- initial state
  double lambda0 = inp["lambda0"], s1 = inp["s1"], rho1 = inp["rho1"];
  Rcpp::NumericVector Rv = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(inp["Rv"]));
  std::vector<double> Wv(N);
  for (int i = 0; i < N; i++) Wv[i] = exp(Rv[i]);

  // the forest: the hypers and options of MakeForest(hyp, Opts(...)) in sim_worker.R.
  // Keeping a copy of every forest (cache_trees) changes no draw and only costs memory.
  Rcpp::List op = Rcpp::clone(opts);
  op["cache_trees"] = false;
  if (!op.containsElementNamed("update_num_tree")) op["update_num_tree"] = false;
  Forest forest(hypers, op);

  CarLD car_ld(N, A.begin(), D.begin(), logdet_fun);
  arma::mat X(Xr.begin(), nt, p, false, true);
  arma::mat Xd(Xdr.begin(), nt, pd, false, true);
  arma::mat Xempty(0, pd);

  // ---- outputs
  Rcpp::NumericVector keep_lam(num_iter), keep_s1(num_iter), keep_rho1(num_iter);
  Rcpp::NumericMatrix keep_W(num_iter, N);
  std::vector<double> acc_R(N, 0.0);
  double acc_rho1 = 0;
  std::vector<double> S_amse_sum((size_t) N * nXP * nT, 0.0), S_aes_sum(4 * nPT, 0.0);
  Rcpp::NumericVector S_aes_draws((size_t) 4 * nPT * n_keep_max, NA_REAL);
  Rcpp::NumericVector S_aesM_draws(use_step3 == 1 ? (size_t) 4 * nPT * n_keep_max : 0, NA_REAL);
  Rcpp::NumericVector lw_keep(n_keep_max);
  std::vector<double> lw_all;
  int n_keep = 0, n_acc = 0;
  double sumw = 0, maxlw = R_NegInf, w_lam = 0, w_rho1 = 0, w_s1 = 0;
  std::vector<double> w_W(N, 0.0);

  // ---- fixed pieces
  std::vector<int> ev_idx;
  for (int k = 0; k < nt; k++) if (delta[k] == 1) ev_idx.push_back(k);
  int nev = ev_idx.size();
  int sum_delta = 0;
  for (int k = 0; k < nt; k++) sum_delta += delta[k];
  std::vector<double> prop(N), Rp(N);
  std::vector<double> ph_buf, ctr_buf, pa((size_t) nPT0 * 4), tpa((size_t) 4 * nPT0), ca(4 * nPT);
  std::vector<double> pq_m, pq_h, py_m((size_t) nt), py_h((size_t) nt);
  Rcpp::NumericMatrix Xq_m = Rcpp::clone(Xq_base), Xq_h = Rcpp::clone(Xq_base);
  Rcpp::NumericMatrix Xy_m = Rcpp::clone(Xy_base), Xy_h = Rcpp::clone(Xy_base);
  Rcpp::NumericMatrix aesM = Rcpp::clone(aesM_base);
  int nq = Xq_base.nrow();
  pq_m.resize(nq); pq_h.resize(nq);
  for (int r = 0; r < nq; r++) Xq_h(r, 1) = M_hat[cty_q[r] - 1];
  for (int k = 0; k < nt; k++) Xy_h(k, 1) = M_hat[cty[k] - 1];
  std::vector<double> Im(nt), Ih(nt);
  size_t nrow_blk = (size_t) nXP * nT0;
  ph_buf.resize(nrow_blk);
  ctr_buf.resize((size_t) nXP * nT);
  std::vector<double> m_cnt(N);

  for (int iter = 1; iter <= num_iter; iter++) {

    // augmentation: the rejected points of the thinned Poisson process
    std::vector<int> q(nt);
    for (int k = 0; k < nt; k++) {
      double lw_ = lambda0 * Wv[cty[k] - 1];
      double rate = lw_ * time[k];
      q[k] = (int) R::rpois(std::max(rate, 1e-8));
    }
    std::vector<int> rep_i;
    for (int k = 0; k < nt; k++) for (int c = 0; c < q[k]; c++) rep_i.push_back(k);
    int nG = rep_i.size();
    std::fill(m_cnt.begin(), m_cnt.end(), 0.0);
    std::vector<int> keepG;
    arma::mat XG;
    if (nG > 0) {
      XG.set_size(nG, pd);
      for (int g = 0; g < nG; g++) {
        int k = rep_i[g];
        double G = R::runif(0, time[k]);
        XG(g, 0) = G / max_time;
        XG(g, 1) = M_hat[cty[k] - 1];
        for (int j = 0; j < p; j++) XG(g, 2 + j) = X(k, j);
      }
      arma::vec bG = forest.do_predict(XG);
      for (int g = 0; g < nG; g++) {
        double u = R::runif(0, 1);
        if (u > R::pnorm(bG(g), 0.0, 1.0, 1, 0)) keepG.push_back(g);
      }
      for (int g : keepG) m_cnt[cty[rep_i[g]] - 1] += 1;
    }
    int naug = keepG.size();
    arma::mat Xall(naug + nev, pd);
    for (int a = 0; a < naug; a++) Xall.row(a) = XG.row(keepG[a]);
    for (int e = 0; e < nev; e++) Xall.row(naug + e) = Xd.row(ev_idx[e]);
    int nall = naug + nev;
    // events are 1, rejected points 0; latent Z, then one back-fitting sweep.
    // As in ifelse(), both rtruncnorm() vectors are drawn in full: first the one
    // truncated to (0, Inf), then the one truncated to (-Inf, 0).
    arma::vec bcur = forest.do_predict(Xall);
    arma::vec Z(nall);
    std::vector<double> yes(nall), no(nall);
    if (nev > 0) for (int k = 0; k < nall; k++) yes[k] = tn_left(0.0, bcur(k), 1.0);
    if (naug > 0) for (int k = 0; k < nall; k++) no[k] = tn_right(0.0, bcur(k), 1.0);
    for (int k = 0; k < nall; k++) Z(k) = (k >= naug) ? yes[k] : no[k];
    // Called from R, SoftBart's forest$do_gibbs() does not save R's random-number
    // state, so the next R draw starts again from the state before the sweep and
    // reuses the sweep's uniforms. To reproduce sim_worker.R the state is saved
    // before the sweep and restored after it.
    PutRNGstate();
    forest.do_gibbs(Xall, Z, Xempty, 1);
    GetRNGstate();

    // lambda0 and sigma1^2 are closed form; rho1 and R need MH
    LDOUBLE swt = 0.0;
    for (int k = 0; k < nt; k++) { double v = Wv[cty[k] - 1] * time[k]; swt += v; }
    double sum_m = 0;
    for (int i = 0; i < N; i++) sum_m += m_cnt[i];
    double shape_l = (1 + sum_delta) + sum_m;
    double rate_l = 1 + (double) swt;
    lambda0 = R::rgamma(shape_l, 1 / rate_l);

    double qf = Rcpp::as<double>(qf_fun(Rv, rho1));   // from R, as in Step 1
    double rate_s = 1 + qf / 2;
    s1 = 1 / R::rgamma(1 + N / 2.0, 1 / rate_s);

    double rp = R::runif(std::max(rb[0], rho1 - 0.05), std::min(rb[1], rho1 + 0.05));
    double u = R::runif(0, 1);
    if (log(u) < car_ld(Rv.begin(), rp, s1) - car_ld(Rv.begin(), rho1, s1)) {
      rho1 = rp; acc_rho1 = acc_rho1 + 1;
    }

    // the W exponent counts events AND augmented points
    for (int c2 = 0; c2 < N; c2++) {
      for (int k = 0; k < N; k++) Rp[k] = Rv[k];
      Rp[c2] = R::rnorm(Rv[c2], 0.2);
      double cnt = d_cnt[c2] + m_cnt[c2];
      double a1 = cnt * Rp[c2], b1 = lambda0 * exp(Rp[c2]), c1 = b1 * sy[c2];
      double llp = a1 - c1;
      double a2 = cnt * Rv[c2], b2 = lambda0 * exp(Rv[c2]), c2b = b2 * sy[c2];
      double llc = a2 - c2b;
      double la = (llp + car_ld(Rp.data(), rho1, s1)) - (llc + car_ld(Rv.begin(), rho1, s1));
      if (R_FINITE(la) && log(R::runif(0, 1)) < la) {
        for (int k = 0; k < N; k++) Rv[k] = Rp[k];
        acc_R[c2] = acc_R[c2] + 1;
      }
    }
    for (int i = 0; i < N; i++) Wv[i] = exp(Rv[i]);
    keep_lam[iter - 1] = lambda0; keep_s1[iter - 1] = s1; keep_rho1[iter - 1] = rho1;
    for (int i = 0; i < N; i++) keep_W(iter - 1, i) = Wv[i];

    if (iter > burn_in && iter % THIN == 0) {
      double lw = 0;
      int mrow = -1;
      // ---- Step 3: the importance weight for this draw
      if (use_step3 == 1) {
        mrow = (int) R_unif_index((double) M_post.nrow());
        for (int r = 0; r < nq; r++) Xq_m(r, 1) = M_post(mrow, cty_q[r] - 1);
        for (int k = 0; k < nt; k++) Xy_m(k, 1) = M_post(mrow, cty[k] - 1);
        predict_mt(forest, Xq_m.begin(), nq, pd, pq_m.data(), nthreads);
        predict_mt(forest, Xq_h.begin(), nq, pd, pq_h.data(), nthreads);
        predict_mt(forest, Xy_m.begin(), nt, pd, py_m.data(), 1);
        predict_mt(forest, Xy_h.begin(), nt, pd, py_h.data(), 1);
        for (int k = 0; k < nt; k++) {
          LDOUBLE am = 0.0, ah = 0.0;
          for (int j = 0; j < NQ; j++) {
            am += R::pnorm(pq_m[k + (size_t) nt * j], 0.0, 1.0, 1, 0);
            ah += R::pnorm(pq_h[k + (size_t) nt * j], 0.0, 1.0, 1, 0);
          }
          am /= NQ; ah /= NQ;
          Im[k] = ((double) am) * time[k];
          Ih[k] = ((double) ah) * time[k];
        }
        LDOUBLE s_a = 0.0, s_b = 0.0;
        for (int k = 0; k < nt; k++) {
          double lPm = log(std::max(R::pnorm(py_m[k], 0.0, 1.0, 1, 0), 1e-300));
          double lPh = log(std::max(R::pnorm(py_h[k], 0.0, 1.0, 1, 0), 1e-300));
          double d1 = lPm - lPh;
          double t1 = (double) delta[k] * d1;
          s_a += t1;
          double d2 = Im[k] - Ih[k];
          double t2 = Wv[cty[k] - 1] * d2;
          s_b += t2;
        }
        double sb_ = (double) s_b;
        double t3 = lambda0 * sb_;
        lw = (double) s_a - t3;
        lw_all.push_back(lw);
      }
      // running max keeps exp() from overflowing
      double w;
      if (lw > maxlw) {
        if (R_FINITE(maxlw)) {
          double f = exp(maxlw - lw);
          for (auto& v : S_amse_sum) v = v * f;
          for (auto& v : S_aes_sum) v = v * f;
          w_lam = w_lam * f;
          for (auto& v : w_W) v = v * f;
          w_rho1 = w_rho1 * f; w_s1 = w_s1 * f; sumw = sumw * f;
        }
        maxlw = lw; w = 1;
      } else w = exp(lw - maxlw);
      n_acc = n_acc + 1; sumw = sumw + w;
      { double t = w * lambda0; w_lam = w_lam + t; }
      for (int i = 0; i < N; i++) { double t = w * Wv[i]; w_W[i] = w_W[i] + t; }
      { double t = w * rho1; w_rho1 = w_rho1 + t; }
      { double t = w * s1; w_s1 = w_s1 + t; }

      for (int i = 0; i < N; i++) {
        Rcpp::NumericMatrix blk = amse_blocks[i];
        predict_mt(forest, blk.begin(), nrow_blk, pd, ph_buf.data(), nthreads);
        for (size_t r = 0; r < nrow_blk; r++) ph_buf[r] = R::pnorm(ph_buf[r], 0.0, 1.0, 1, 0);
        cumtrap(ph_buf.data(), nXP, nT0, dt, ctr_buf.data());
        double sc = -lambda0 * Wv[i];
        for (int r = 0; r < nXP; r++) for (int k = 0; k < nT; k++) {
          double e = sc * ctr_buf[r + (size_t) nXP * k];
          double v = clip01(exp(e));
          double t = w * v;
          size_t id = i + (size_t) N * (r + (size_t) nXP * k);
          S_amse_sum[id] = S_amse_sum[id] + t;
        }
      }
      predict_mt(forest, aes_block.begin(), aes_block.nrow(), pd, pa.data(), 1);
      for (int j = 0; j < nPT0; j++) for (int c = 0; c < 4; c++)
        tpa[c + 4 * j] = R::pnorm(pa[j + (size_t) nPT0 * c], 0.0, 1.0, 1, 0);
      cumtrap(tpa.data(), 4, nPT0, dtp, ca.data());
      double sc1 = -lambda0 * Wv[0];
      n_keep = n_keep + 1;
      for (int idx = 0; idx < 4 * nPT; idx++) {
        double e = sc1 * ca[idx];
        double v = clip01(exp(e));
        double t = w * v;
        S_aes_sum[idx] = S_aes_sum[idx] + t;
        S_aes_draws[(size_t) 4 * nPT * (n_keep - 1) + idx] = v;
      }
      lw_keep[n_keep - 1] = lw;

      if (use_step3 == 1) {
        for (int r = 0; r < aesM.nrow(); r++) aesM(r, 1) = M_post(mrow, 0);
        predict_mt(forest, aesM.begin(), aesM.nrow(), pd, pa.data(), 1);
        for (int j = 0; j < nPT0; j++) for (int c = 0; c < 4; c++)
          tpa[c + 4 * j] = R::pnorm(pa[j + (size_t) nPT0 * c], 0.0, 1.0, 1, 0);
        cumtrap(tpa.data(), 4, nPT0, dtp, ca.data());
        for (int idx = 0; idx < 4 * nPT; idx++) {
          double e = sc1 * ca[idx];
          S_aesM_draws[(size_t) 4 * nPT * (n_keep - 1) + idx] = clip01(exp(e));
        }
      }
    }
    if (iter % 1000 == 0) {
      log_iter(iter, lambda0, s1, rho1, Rcpp::wrap(acc_R), acc_rho1);
    }
    if (iter % 100 == 0) Rcpp::checkUserInterrupt();
  }

  Rcpp::NumericVector S_amse(S_amse_sum.begin(), S_amse_sum.end());
  S_amse.attr("dim") = Rcpp::IntegerVector::create(N, nXP, nT);
  Rcpp::NumericMatrix S_aes(4, nPT);
  std::copy(S_aes_sum.begin(), S_aes_sum.end(), S_aes.begin());
  S_aes_draws.attr("dim") = Rcpp::IntegerVector::create(4, nPT, n_keep_max);
  Rcpp::RObject aesM_out = R_NilValue;
  if (use_step3 == 1) {
    S_aesM_draws.attr("dim") = Rcpp::IntegerVector::create(4, nPT, n_keep_max);
    aesM_out = S_aesM_draws;
  }
  Rcpp::RObject lw_all_out = R_NilValue;
  if (!lw_all.empty()) lw_all_out = Rcpp::NumericVector(lw_all.begin(), lw_all.end());

  Rcpp::List out;
  out["keep_lam"] = keep_lam; out["keep_s1"] = keep_s1; out["keep_rho1"] = keep_rho1;
  out["keep_W"] = keep_W; out["acc_R"] = Rcpp::wrap(acc_R); out["acc_rho1"] = acc_rho1;
  out["S_amse_sum"] = S_amse; out["S_aes_sum"] = S_aes; out["sumw"] = sumw;
  out["maxlw"] = maxlw; out["lw_all"] = lw_all_out; out["w_lam"] = w_lam;
  out["w_W"] = Rcpp::wrap(w_W); out["w_rho1"] = w_rho1; out["w_s1"] = w_s1;
  out["n_acc"] = (double) n_acc; out["S_aes_draws"] = S_aes_draws; out["S_aesM_draws"] = aesM_out;
  out["lw_keep"] = lw_keep; out["n_keep"] = (double) n_keep; out["lambda0"] = lambda0;
  out["Rv"] = Rv; out["Wv"] = Rcpp::wrap(Wv); out["s1"] = s1; out["rho1"] = rho1;
  return out;
}
