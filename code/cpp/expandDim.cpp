// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../CharlieLib/cpp/Charlie/VecPool.hpp"
#include "../CharlieLib/cpp/Charlie/ThreadPool.hpp"
#include "../CharlieLib/cpp/Charlie/ProgressBar.hpp"
#define vec std::vector


void eventLossInAggCounties(
    int *ind, double *loss, int size,
    int *countyInd, int countyIndSize,
    double *rst) // rst is of size countyIndSize.
{
  // std::fill(rst, rst + countyIndSize, 0.0);
  if constexpr (true)
  {
    constexpr const int i = 0;
    auto it = std::lower_bound(ind, ind + size, countyInd[i]);
    if (it != ind + size and *it == countyInd[i]) // Found.
      rst[i] = loss[it - ind];
    else rst[i] = 0;
  }
  for (int i = 1; i < countyIndSize; ++i)
  {
    auto it = std::lower_bound(ind, ind + size, countyInd[i]);
    if (it != ind + size and *it == countyInd[i])
      rst[i] = rst[i - 1] + loss[it - ind];
    else rst[i] = rst[i - 1];
  }
}


// Given an event, return the extended event.
std::pair<vec<int>, vec<float>> f(
    int *ind, double *loss, int size, 
    int* nmat, int nrow, int ncol,
    int *levelInterest, int Nlevel, // e.g. 8, 16, 32, 64
    bool takeAvg,
    int dimStart, // First index that can be used.
    Charlie::VecPool &vp)
{
  auto cumuLoss = vp.lend<double> (nrow * (int64_t)ncol);
  for (int i = 0; i < ncol; ++i)
  {
    std::size_t offset = nrow * (std::size_t)i;
    eventLossInAggCounties(
      ind, loss, size, nmat + offset, nrow, cumuLoss.data() + offset);
  }
  auto Vind = vp.lend<int> (Nlevel, ncol);
  for (auto &x: Vind) x.resize(0);
  auto Vloss = vp.lend<float> (Nlevel, ncol);
  for (auto &x: Vloss) x.resize(0);
  
  
  for (int i = 0; i < ncol; ++i)
  {
    for (int j = 0; j < Nlevel; ++j)
    {
      int l = levelInterest[j] - 1;
      // if ( cumuLoss[i][l] > 1e-10 )
      double lss = cumuLoss[i * nrow + l];
      if ( lss > 1e-10 )
      {
        // Vind[j].emplace_back(  i + dimStartV[j]  );
        Vind[j].emplace_back(  i + dimStart + ncol * j  );
        // Vloss[j].emplace_back( cumuLoss[i][l]  );
        // double lss = cumuLoss[i * nrow + l];
        lss = takeAvg ? lss / (l + 1) : lss;
        Vloss[j].emplace_back( lss);
      }
    }
  } 
  
  
  std::pair<vec<int>, vec<float>> rst;
  int Ndim = size;
  for (auto &x: Vind) Ndim += x.size();
  rst.first.reserve(Ndim);
  rst.second.reserve(Ndim);
  rst.first.assign(ind, ind + size);
  rst.second.assign(loss, loss + size);
  for (int i = 0; i < Nlevel; ++i)
  { 
    for (int j = 0, jend = Vind[i].size(); j < jend; ++j)
    { 
      rst.first.emplace_back(Vind[i][j]);
      rst.second.emplace_back(Vloss[i][j]);
    }
  } 
  
  
  vp.recall(Vloss);
  vp.recall(Vind);
  vp.recall(cumuLoss);
  return rst;
} 




// NNs = e.g. c(8, 16, 32, 64). nmat should have at least 64 rows.
// [[Rcpp::export]]
List makeExpandedEvents(List eventLosses, IntegerMatrix nmat,
                        int dimStart,
                        IntegerVector NNs, bool takeAvg = true,
                        int maxCore = 10000, bool verbose = true) 
{ 
  struct event { int size, *ind; double *loss; };
  vec<event> E(eventLosses.size());
  // int dimStart = 0;
  for (int i = 0, iend = E.size(); i < iend; ++i)
  { 
    List l = eventLosses[i];
    IntegerVector ind = l[0];
    NumericVector loss = l[1];
    E[i].ind = &ind[0];
    E[i].loss = &loss[0];
    E[i].size = ind.size();
    // for (auto &x: ind) dimStart = std::max(dimStart, x);
  } 
  // dimStart += 1;
  
  
  vec<std::pair<vec<int>, vec<float> > > rst(E.size());
  Charlie::ThreadPool cp(std::move(maxCore));
  vec<Charlie::VecPool> vps(maxCore);
  Charlie::ProgressBar<100> pb(E.size());
  if (verbose) Rcout << "Progress %: ";
  auto g = [&](std::size_t i, std::size_t t)->bool
  { 
    if (verbose and t == 0)
    {
      auto p = pb(i);
      if (p != -1) Rcout << p << " ";
    }
    auto &vp = vps[t];
    auto vecPair = f(E[i].ind, E[i].loss, E[i].size, &nmat[0], nmat.nrow(), 
                     nmat.ncol(), &NNs[0], NNs.size(), takeAvg, dimStart, vp);
    rst[i] = std::move(vecPair);
    return false;
  };
  cp.parFor(0, E.size(), g, E.size() / (maxCore * maxCore * maxCore) + 1);
  
  
  List result(E.size());
  for (int i = 0, iend = E.size(); i < iend; ++i)
  {
    IntegerVector ind(rst[i].first.begin(), rst[i].first.end());
    vec<int>().swap(rst[i].first);
    NumericVector loss(rst[i].second.begin(), rst[i].second.end());
    vec<float>().swap(rst[i].second);
    result[i] = List::create(ind, loss);
  }
  return result;
} 







#undef vec




