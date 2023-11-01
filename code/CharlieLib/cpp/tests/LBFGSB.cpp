// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie.hpp"
#define vec std::vector


struct Rosenbrock
{
  double operator()(double *x, double *&grad, int dim)
  {
    // grad.resize(dim);
    double fx = 0.0;
    for(int i = 0; i < dim; i += 2)
    {
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      grad[i + 1] = 20 * t2;
      grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
      fx += t1 * t1 + t2 * t2;
    } 
    return fx;
  }
};


struct RosenbrockNoGrad
{
  double operator()(double *x, double *&grad, int dim)
  {
    // grad.resize(dim);
    double fx = 0.0;
    for(int i = 0; i < dim; i += 2)
    { 
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      fx += t1 * t1 + t2 * t2;
    }
    grad = nullptr;
    return fx;
  } 
};




// [[Rcpp::export]]
List miniRosenbrock(NumericVector x, int maxit = 100, bool giveGrad = false)
{
  if (x.size() % 2 != 0) stop("x.size() is not even!");
  NumericVector rst(x.begin(), x.end());
  Charlie::LBFGSB lbfgsb;
  if (giveGrad)
  {
    Rosenbrock f;
    auto fval = lbfgsb(&rst[0], x.size(), f, maxit);
    return List::create(Named("param") = rst, 
                        Named("fval") = fval.first,
                        Named("Niter") = fval.second);  
  }
  RosenbrockNoGrad f;
  auto fval = lbfgsb(&rst[0], x.size(), f, maxit);
  return List::create(Named("param") = rst, 
                      Named("fval") = fval.first,
                      Named("Niter") = fval.second);
}









