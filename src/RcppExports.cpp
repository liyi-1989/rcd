// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ccorecpp
double ccorecpp(NumericVector u, NumericVector v, NumericVector bw);
RcppExport SEXP rcd_ccorecpp(SEXP uSEXP, SEXP vSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    __result = Rcpp::wrap(ccorecpp(u, v, bw));
    return __result;
END_RCPP
}
// ccorecppabs
double ccorecppabs(NumericVector u, NumericVector v, NumericVector bw);
RcppExport SEXP rcd_ccorecppabs(SEXP uSEXP, SEXP vSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    __result = Rcpp::wrap(ccorecppabs(u, v, bw));
    return __result;
END_RCPP
}
// rcpp_distance
NumericMatrix rcpp_distance(NumericMatrix mat);
RcppExport SEXP rcd_rcpp_distance(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    __result = Rcpp::wrap(rcpp_distance(mat));
    return __result;
END_RCPP
}
// rcpp_parallel_distance
NumericMatrix rcpp_parallel_distance(NumericMatrix mat);
RcppExport SEXP rcd_rcpp_parallel_distance(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    __result = Rcpp::wrap(rcpp_parallel_distance(mat));
    return __result;
END_RCPP
}
// Kcpp
double Kcpp(double u);
RcppExport SEXP rcd_Kcpp(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    __result = Rcpp::wrap(Kcpp(u));
    return __result;
END_RCPP
}
// PKcpp
double PKcpp(NumericVector x, NumericVector xi, NumericVector h);
RcppExport SEXP rcd_PKcpp(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    __result = Rcpp::wrap(PKcpp(x, xi, h));
    return __result;
END_RCPP
}
// kdendcpp
double kdendcpp(NumericVector x, NumericMatrix X, NumericVector h);
RcppExport SEXP rcd_kdendcpp(SEXP xSEXP, SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    __result = Rcpp::wrap(kdendcpp(x, X, h));
    return __result;
END_RCPP
}
// kdendveccpp
NumericVector kdendveccpp(NumericMatrix X, NumericVector h);
RcppExport SEXP rcd_kdendveccpp(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    __result = Rcpp::wrap(kdendveccpp(X, h));
    return __result;
END_RCPP
}
