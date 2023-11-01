// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie.hpp"
#define vec std::vector // Always undefined it at file end.


// [[Rcpp::export]]
IntegerVector test(int seed, int someInt) 
{
  // Charlie::SAAmemPool sap(300);
  Charlie::MiniPCG rng(seed);
  std::uniform_real_distribution<double> U(0, 1);
  vec<char> rst;
  auto push = [&](char* x, char* xend)->void
  {
    for (; x < xend; ++x) 
      rst.emplace_back(*x);
  };
  
  
  Charlie::SAAmemPool mempool;
  auto X = Charlie::mvecMake<char> (mempool, rng() % someInt + someInt, 0);
  
  
  X.resize(0);
  for (int i = 0; i < someInt; ++i)
  {
    double u = U(rng);
    if (u < 1.0 / 3 and X.size() != 0)
    {
      int k = rng() % X.size();
      auto& x = X[k];
      push(&*x.begin(), &*x.end());
      // vec<char>(0).swap(x);
      // Charlie::mvecMake<char> (0, mempool).swap(x);
      Charlie::mvecMake<char> (mempool, 0).swap(x);
    }
    else if (u < 2.0 / 3)
    {
      // vec<char> tmp(rng() % someInt);
      // auto tmp = Charlie::mvecMake<char>(rng() % someInt, mempool);
      auto tmp = Charlie::mvecMake<char>(mempool, rng() % someInt);
      for (auto &x: tmp) x = char(rng() % 256);
      X.emplace_back(std::move(tmp));
    }
    else if (X.size() != 0)
    {
      int k = rng() % X.size();
      auto &x = X[k];
      for (int t = 0; t < someInt; ++t)
        x.emplace_back(char(rng() % 256));
    }
  }
  
  
  for (auto& x: X) push(&*x.begin(), &*x.end());
  
  
  return IntegerVector(rst.begin(), rst.end());
}



// [[Rcpp::export]]
IntegerVector testBenchmark(int seed, int someInt) 
{
  // Charlie::SAAmemPool sap(300);
  Charlie::MiniPCG rng(seed);
  std::uniform_real_distribution<double> U(0, 1);
  vec<char> rst;
  auto push = [&](char* x, char* xend)->void
  {
    for (; x < xend; ++x) 
      rst.emplace_back(*x);
  };
  
  
  vec<vec<char>> X(rng() % someInt + someInt);
  X.resize(0);
  for (int i = 0; i < someInt; ++i)
  {
    double u = U(rng);
    if (u < 1.0 / 3 and X.size() != 0)
    {
      int k = rng() % X.size();
      auto& x = X[k];
      push(&*x.begin(), &*x.end());
      vec<char>(0).swap(x);
    }
    else if (u < 2.0 / 3)
    {
      vec<char> tmp(rng() % someInt);
      for (auto &x: tmp) x = char(rng() % 256);
      X.emplace_back(std::move(tmp));
    }
    else if (X.size() != 0)
    {
      int k = rng() % X.size();
      auto &x = X[k];
      for (int t = 0; t < someInt; ++t)
        x.emplace_back(char(rng() % 256));
    }
  }
  
  
  for (auto& x: X) push(&*x.begin(), &*x.end());
  
  
  return IntegerVector(rst.begin(), rst.end());
}






#undef vec




