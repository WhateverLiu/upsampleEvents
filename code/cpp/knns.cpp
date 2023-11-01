// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "hpp/distance.hpp"
#include "../CharlieLib/cpp/Charlie/ThreadPool.hpp"
#include "../CharlieLib/cpp/Charlie/ProgressBar.hpp"
#define vec std::vector


template <typename ing, typename num, typename Dtype, bool useWeight>
void knnsCore(
    Event<ing, num> *event, const ing Nevent,
    const ing NcoreEvent,
    ing *candiMat, const ing NcandiEach, // NcandiEach = nrow(candiMat)
    ing *rst,
    const int keep, 
    int maxCore, 
    const bool verbose,
    num *w, double p)
{
  vec<vec<std::pair<double, ing> > > pV;
  Charlie::ProgressBar<100> pb(Nevent - NcoreEvent);
  if (verbose) Rcout << "Progress %: ";
  
  
  // const ing NirreEvent = Nevent - NcoreEvent; // Number of events to be replaced.
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    if (verbose and t == 0)
    {
      auto p = pb(i - NcoreEvent);
      if (p != -1) Rcout << p << " ";
    }
    auto &e = event[i];
    const ing *cans = (i - NcoreEvent) * std::size_t(NcandiEach) + candiMat;
    auto &cntr = pV[t];
    cntr.resize(NcandiEach);
    Distance<ing, num, Dtype, useWeight> dst;
    if constexpr (std::is_same<Dtype, LpUnit>::value) dst.dist.p = p;
    for (int k = 0; k < NcandiEach; ++k)
    {
      cntr[k].first = dst(e, event[cans[k]], w);
      cntr[k].second = event[cans[k]].id;
    }
    std::partial_sort(cntr.begin(), cntr.begin() + keep, cntr.end());
    for (int k = 0; k < keep; ++k) 
      rst[(i - NcoreEvent) * keep + k] = cntr[k].second;
    return false;
  };
  Charlie::ThreadPool cp(std::move(maxCore));
  pV.resize(maxCore);
  cp.parFor(NcoreEvent, Nevent, f);
}




template <typename ing, typename num, typename Dtype>
void knns(
    Event<ing, num> *event, const ing Nevent,
    const ing NcoreEvent,
    ing *candiMat, const ing NcandiEach, // NcandiEach = nrow(candiMat)
    ing *rst,
    const int keep, 
    int maxCore, 
    const bool verbose,
    num *w, double p)
{
  if (w == nullptr)
    knnsCore<ing, num, Dtype, false>(
      event, Nevent,
      NcoreEvent,
      candiMat, NcandiEach, // NcandiEach = nrow(candiMat)
      rst,
      keep, 
      maxCore, 
      verbose,
      w, p);
  else
    knnsCore<ing, num, Dtype, true>(
        event, Nevent,
        NcoreEvent,
        candiMat, NcandiEach, // NcandiEach = nrow(candiMat)
        rst,
        keep, 
        maxCore, 
        verbose,
        w, p);
}
  



// eventLosses.size() - candidates.size() is the number of the events in the
//   first 10K years.
// candidates are candidate events for the last 90K catalog.
// [[Rcpp::export]]
IntegerMatrix knns(List eventLosses, IntegerMatrix candidateMat, 
                   NumericVector regionW = NumericVector(0),
                   int keep = 100, 
                   String distance = "Euclidean", 
                   double LPp = 2, // If LP measure is used, LPp is the power.
                   int maxCore = 10000, 
                   bool verbose = true) 
{
  if (std::abs(LPp - 2) < 1e-10 and distance == "LP") distance = "Euclidean";
  else if (std::abs(LPp - 1) < 1e-10 and distance == "LP") distance = "L1";
  else if (std::abs(LPp - 0.5) < 1e-10 and distance == "LP") distance = "L-half";
  else if (std::abs(LPp - 0.25) < 1e-10 and distance == "LP") distance = "quarterRoot";
  else if (std::abs(LPp - 1.0 / 3) < 1e-10 and distance == "LP") distance = "cubitRoot";
  
  
  const int NcoreEvent = eventLosses.size() - candidateMat.ncol();
  if (candidateMat.nrow() < keep) stop("keep is too large.");
  
  
  vec<Event<int, double> > events(eventLosses.size());
  std::size_t Nslots = 0;
  for (int i = 0, iend = events.size(); i < iend; ++i)
  {
    List l = eventLosses[i];
    IntegerVector ind = l[0];
    NumericVector val = l[1];
    events[i].reset(&ind[0], &val[0], i, ind.size());
    Nslots += ind.size();
  }
  
  
  vec<double> binarySlots(Nslots);
  if (distance == "binary")
  {
    std::size_t u = 0;
    for (int i = 0, iend = events.size(); i < iend; ++i)
    {
      double *v = binarySlots.data() + u;
      for (int k = 0, kend = events[i].size; k < kend; ++k)
      {
        binarySlots[u] = events[i].val[k] > 1e-10;
        u += 1;
      }
      events[i].val = v;
    }
    distance = "Euclidean";
  }
  
  
  double *w = nullptr;
  if (regionW.size() != 0) w = &regionW[0];
  
  
  IntegerMatrix rst(keep, candidateMat.ncol());
  
  
  if (distance == "Euclidean")
  {
    for (auto &x: events)
    {
      if (w != nullptr)
      {
        x.sumOfSquaredVal = 0;
        for (int i = 0; i < x.size; ++i)
          x.sumOfSquaredVal += x.val[i] * x.val[i] * w[x.ind[i]];  
      }
      else x.sumOfSquaredVal = std::inner_product(x.val, x.val + x.size, x.val, 0.0); 
    }
    knns<int, double, EucDunit>(
        events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
        &rst[0], keep, maxCore, verbose, w, LPp); 
  }
  else if (distance == "L1") knns<int, double, L1unit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "L-half") knns<int, double, LhalfUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "symCrossEntropy") knns<int, double, EntropyUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "dotProduct") knns<int, double, DotProdUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "cubitRoot") knns<int, double, LcubicRootUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "quarterRoot") knns<int, double, LquarterUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else if (distance == "LP") knns<int, double, LpUnit>(
    events.data(), events.size(), NcoreEvent, &candidateMat[0], candidateMat.nrow(), 
    &rst[0], keep, maxCore, verbose, w, LPp);
  else stop("Distance measure not implemented.");
  
  
  return rst;
}




#undef vec









