#include "Event.hpp"


struct DotProdUnit
{
  double operator()(double x, double y, double w) {  return -x * y * w; }
  constexpr double operator()(double x, double w) {  return 0; } 
  double output(double S) { return S; }
};


// Will be different from DotProdUnit in the end.
struct EucDunit: public DotProdUnit 
{
  double output(double S) { return std::sqrt(S); }
};


struct L1unit
{
  double operator()(double x, double y, double w) { return std::abs(x - y) * w; }
  double operator()(double x, double w) { return std::abs(x) * w; }
  double output(double S) { return S; }
};


struct LhalfUnit
{
  double operator()(double x, double y, double w) { return std::sqrt(std::abs(x - y)) * w; }
  double operator()(double x, double w) { return std::sqrt(std::abs(x)) * w; }
  double output(double S) { return S * S; }
}; 


struct LquarterUnit
{
  double operator()(double x, double y, double w) { return std::sqrt(std::sqrt(std::abs(x - y))) * w; }
  double operator()(double x, double w) { return std::sqrt(std::sqrt(std::abs(x))) * w; }
  double output(double S) { return (S * S) * (S * S); }
};  


struct LcubicRootUnit
{
  double operator()(double x, double y, double w) { return std::cbrt(std::abs(x - y)) * w; }
  double operator()(double x, double w) { return std::cbrt(std::abs(x)) * w; }
  double output(double S) { return S * S * S; }
};   


struct EntropyUnit
{
  double operator()(double x, double y, double w) 
  { 
    return -(x * std::log(y) + y * std::log(x)) * w;
  }
  double operator()(double x, double w)
  { 
    constexpr const double smallVal = 1e-300;
    return -(x * std::log(smallVal) + smallVal * std::log(x)) * w;
  }
  double output(double S) { return S; }
};


struct LpUnit
{
  double p;
  double operator()(double x, double y, double w) { return std::pow(std::abs(x - y), p) * w; }
  double operator()(double x, double w) { return std::pow(std::abs(x), p) * w; }
  double output(double S) { return std::pow(S, 1.0 / p); }
};


template<typename ing, typename num, typename Dunit, bool useWeight>
struct Distance
{
  Dunit dist;
  double operator()(Event<ing, num> &x, Event<ing, num> &y, num *w)
  {
    ing i = 0, j = 0;
    double rst = 0;
    while (i < x.size and j < y.size)
    {
      if (x.ind[i] == y.ind[j])
      {
        if (useWeight) rst += dist(x.val[i], y.val[j], w[x.ind[i]]);
        else rst += dist(x.val[i], y.val[j], 1.0);
        i += 1; j += 1;
      }
      else if (x.ind[i] < y.ind[j])
      {
        if (useWeight) rst += dist(x.val[i], w[x.ind[i]]);
        else rst += dist(x.val[i], 1.0);
        i += 1;
      }
      else
      {
        if (useWeight) rst += dist(y.val[j], w[y.ind[j]]);
        else rst += dist(y.val[j], 1.0);
        j += 1;
      }
    }
    for (; i < x.size; ++i)
    {
      if (useWeight) rst += dist(x.val[i], w[x.ind[i]]); 
      else rst += dist(x.val[i], 1.0); 
    }
    for (; j < y.size; ++j)
    {
      if (useWeight) rst += dist(y.val[j], w[y.ind[j]]); 
      else rst += dist(y.val[j], 1.0); 
    }
    
    
    if ( std::is_same<Dunit, EucDunit>::value )
      rst = std::max(0.0, 2 * rst + x.sumOfSquaredVal + y.sumOfSquaredVal);
    return dist.output(rst);
  }
};



















