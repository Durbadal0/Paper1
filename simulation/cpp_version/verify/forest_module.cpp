// forest_module.cpp -- exposes the engine's copy of the SoftBart forest to R, for
// verify/check_forest_identity.R
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../code/cpp/softbart/soft_bart.h"
#include "../code/cpp/softbart/soft_bart_impl.h"
RCPP_MODULE(engine_forest) {
  Rcpp::class_<Forest>("EngineForest")
    .constructor<Rcpp::List, Rcpp::List>()
    .method("do_gibbs", &Forest::do_gibbs)
    .method("do_predict", &Forest::do_predict)
    .method("get_counts", &Forest::get_counts);
}
