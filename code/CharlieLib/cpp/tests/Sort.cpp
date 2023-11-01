// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie.hpp"


// [[Rcpp::export]]
NumericVector Sort(NumericVector x, int maxCore = 1000) 
{
  NumericVector rst(x.begin(), x.end());
  Charlie::ThreadPool tp(std::move(maxCore));
  Charlie::VecPool vp;
  Charlie::Sort()(rst.begin(), rst.end(), &tp, &vp);
  return rst;
}


/***R

# Rcpp::sourceCpp("cpp/tests/Sort.cpp")
getwd()
tmp = runif(1e8)  
system.time({tmp2 = Sort(tmp)})
system.time({tmp3 = sort(tmp)})
range(tmp3 - tmp2)



*/














