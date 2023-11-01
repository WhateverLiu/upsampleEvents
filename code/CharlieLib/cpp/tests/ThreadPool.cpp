// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie/ThreadPool.hpp"


// [[Rcpp::export]]
int paraSummation(IntegerVector x, int maxCore = 15, int grainSize = 100)
{
  Charlie::ThreadPool cp(std::move(maxCore));
  IntegerVector S(maxCore);
  cp.parFor(0, x.size(), [&](std::size_t i, std::size_t t)->bool
  {
    S[t] += (x[i] % 31 + x[i] % 131 + x[i] % 73 + x[i] % 37 + x[i] % 41) % 7;
    return false; // Return true indicates early return.
  }, grainSize,
  [](std::size_t t)->bool{ return false; },
  [](std::size_t t)->bool{ return false; });
  return std::accumulate(S.begin(), S.end(), 0);
}


/***R
tmp = sample(1000L, 3e7, replace= T)
system.time({cppRst = paraSummation(tmp, maxCore = 100, grainSize = 100)})
system.time({truth = sum((tmp %% 31L + tmp %% 131L + tmp %% 73L +
                            tmp %% 37L + tmp %% 41L) %% 7L)})
cppRst - truth
*/


// [[Rcpp::export]]
double paraSummationFloat(NumericVector x, int maxCore = 15, int grainSize = 100)
{
  Charlie::ThreadPool cp(std::move(maxCore));
  NumericVector S(maxCore);
  cp.parFor(0, x.size(), [&](std::size_t i, std::size_t t)->bool
  {
    S[t] += (x[i] / 31 + x[i] / 131 + x[i] / 73 + x[i] / 37 + x[i] / 41) / 7;
    return false; // Return true indicates early return.
  }, grainSize,
  [](std::size_t t)->bool{ return false; },
  [](std::size_t t)->bool{ return false; });
  return std::accumulate(S.begin(), S.end(), 0.0);
}


/***R
tmp2 = runif(1e8); paraSummationFloat(tmp2, maxCore = 1, grainSize = 100) /
  paraSummationFloat(tmp2, maxCore = 100, grainSize = 100) - 1
*/




