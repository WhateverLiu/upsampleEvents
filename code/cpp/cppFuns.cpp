// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../CharlieLib/cpp/Charlie/Sort.hpp"
#include "../CharlieLib/cpp/Charlie/tiktok.hpp"
#include "../CharlieLib/cpp/Charlie/VecPool.hpp"
#include "../CharlieLib/cpp/Charlie/ThreadPool.hpp"
#define vec std::vector


// [[Rcpp::export]]
double cod(List newEPs, List trueEPs, int catalogSize = 100000, 
           bool onlyAccountForNonzero = false) 
{

  // Charlie::ThreadPool cp(std::move(1));
  // vec<double> error(maxCore, 0);
  double error = 0;
  // vec<double> Sum(maxCore, 0);
  double Sum = 0;
  // vec<double> SumOfSquare(maxCore, 0);
  double SumOfSquare = 0;
  // auto f = [&](std::size_t i, std::size_t t)->bool
  // {
  // NumericMatrix M(3, newEPs.size());
  int Nnonzero = 0;
  for (int i = 0, iend = newEPs.size(); i < iend; ++i)
  {
    NumericVector nep = newEPs[i], tep = trueEPs[i];
    int size = std::min<int> (nep.size(), tep.size());
    Nnonzero += std::max<int> (nep.size(), tep.size());
    double e = 0;
    for (int k = size - 1; k >= 0; --k)
    {
      double tmp = nep[k] - tep[k];
      e += tmp * tmp;
    }
    for (int k = nep.size() - 1; k >= size; --k) e += nep[k] * nep[k];
    for (int k = tep.size() - 1; k >= size; --k) e += tep[k] * tep[k];
    error += e;
    for (int k = tep.size() - 1; k >= 0; --k)
    {
      Sum += tep[k];
      SumOfSquare += tep[k] * tep[k];
    }
    // return false;
  // };
  }
  
  // cp.parFor(0, newEPs.size(), f, 1);
  double Nobs = double(newEPs.size()) * catalogSize;
  if (onlyAccountForNonzero) Nobs = Nnonzero;
  // double EX = std::accumulate(Sum.begin(), Sum.end(), 0.0) / Nobs;
  double EX = Sum / Nobs;
  
  // double EX2 = std::accumulate(SumOfSquare.begin(), SumOfSquare.end(), 0.0) / Nobs;
  double EX2 = SumOfSquare / Nobs;
  double denominator = EX2 - EX * EX;
  // double numerator = std::accumulate(error.begin(), error.end(), 0.0) / Nobs;
  double numerator = error / Nobs;
  return 1.0 - numerator / denominator;
  
  
}



// [v, vend) have been sorted.
void knn1d(std::pair<double, int> *v, 
           std::pair<double, int> *vend, 
           double val, const bool useRelaErr, 
           int *rst, int K)
{
  auto err = [useRelaErr](double x, double y)->double
  {
    if (useRelaErr) return std::abs((x - y) / (x + y));
    return std::abs(x - y);
  };
  
  
  auto it = std::lower_bound(v, vend, val, [](
    const std::pair<double, int> &x, double y)->bool
  {
    return x.first < y;
  });
  
  
  int vsize = vend - v;
  
  
  if (it >= vend)
  {
    int j = vsize - 1;
    for (int i = 0; i < K; ++i, --j) rst[i] = v[j].second;
  }
  else if (it <= v)
  {
    for (int i = 0; i < K; ++i) rst[i] = v[i].second;
  }
  else
  {
    int i = it - v - 1, j = i + 1, k = 0;
    auto ierr = err(v[i].first, val), jerr = err(v[j].first, val);
    while (k < K)
    {
      if (ierr < jerr)
      {
        rst[k] = v[i].second;
        i -= 1;
        k += 1;
        if (i < 0) break;
        ierr = err(v[i].first, val);
      }
      else
      {
        rst[k] = v[j].second;
        j += 1;
        k += 1;
        if (j >= vsize) break;
        jerr = err(v[j].first, val);
      }
    }
    for (; k < K and i >= 0;    --i, ++k) rst[k] = v[i].second;
    for (; k < K and j < vsize; ++j, ++k) rst[k] = v[j].second;
  }
}


// Produce a matrix where each element is 0-based event index.
// [[Rcpp::export]]
IntegerMatrix findCandidates(NumericVector poolEventLosses, 
                        NumericVector lossesOfEventsToBeReplaced,
                        int K = 1000, bool relativeDiff = true, 
                        int maxCore = 10000)
{
  if (K > poolEventLosses.size()) stop("K is too high.");
  // List rst(lossesOfEventsToBeReplaced.size());
  IntegerMatrix rst(K, lossesOfEventsToBeReplaced.size());
  // vec<int*> rstPtr(rst.size());
  // for (int i = 0, iend = rst.size(); i < iend; ++i)
  // {
  //   IntegerVector v(K);
  //   // rst[i] = v;
  //   // rstPtr[i] = &v[0];
  // }
  Charlie::ThreadPool cp(std::move(maxCore));
  vec<std::pair<double, int>> eventLoss(poolEventLosses.size());
  for (int i = 0, iend = eventLoss.size(); i < iend; ++i)
  {
    eventLoss[i].first = poolEventLosses[i];
    eventLoss[i].second = i;
  }
  std::sort(eventLoss.begin(), eventLoss.end());
  const bool useRela = relativeDiff;
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    knn1d(&*eventLoss.begin(), &*eventLoss.end(),
          lossesOfEventsToBeReplaced[i], useRela, 
          // rstPtr[i],
          &rst[0] + K * (std::size_t)i,
                K);
    return false;
  };
  cp.parFor(0, lossesOfEventsToBeReplaced.size(), f, 
            lossesOfEventsToBeReplaced.size() / (maxCore * maxCore * maxCore) + 1);
  return rst;
}




// This shows that hand-crafted makeEPs is not worth it. The major overhead
//   is allocating large list of large vectors in R.
// eventYear can be one-based or zero-based.
// [[Rcpp::export]]
List makeEPs(List eventLosses, IntegerVector eventYear, int maxCore = 10000)
{
  struct Et { int size, *ind; double *val; };
  vec<Et> events(eventLosses.size());
  int64_t unlistSize = 0;
  
  
  Charlie::tiktok< std::chrono::microseconds> timer;
  timer.tik();
  for (int i = 0, iend = events.size(); i < iend; ++i)
  {
    List l = eventLosses[i];
    IntegerVector ind = l[0];
    NumericVector val = l[1];
    events[i].ind = &ind[0];
    events[i].val = &val[0];
    events[i].size = ind.size();
    unlistSize += ind.size();
  }
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  struct A { int region, year; double loss; };
  vec<A> v(unlistSize);
  int Nregion = 0;
  
  
  if constexpr (true)
  {
    std::size_t k = 0;
    for (int i = 0, iend = events.size(); i < iend; ++i)
    {
      for (int j = 0; j < events[i].size; ++j)
      {
        v[k].region = events[i].ind[j];
        Nregion = std::max(Nregion, v[k].region);
        v[k].year = eventYear[i];
        v[k].loss = events[i].val[j];
        k += 1;
      }
    }
    Nregion += 1;
  }
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  Charlie::ThreadPool cp(std::move(maxCore));
  Charlie::Sort()(v.begin(), v.end(), [](const A&x, const A&y)->bool
  {
    if (x.region < y.region) return true;
    if (x.region == y.region) return x.year < y.year;
    return false;
    // return (uint64_t(x.region) << 32) + x.year < (uint64_t(y.region) << 32) + y.year;
  }, &cp);
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  // Get how many losses each region has.
  vec<int> regionNlosses(Nregion);
  regionNlosses.resize(1);
  if constexpr (true)
  {
    regionNlosses.back() = 0;
    // std::size_t prior = 0;
    std::size_t k = 1, end = v.size();
    while (true)
    {
      if (k >= end)
      {
        // regionNlosses.emplace_back(k - prior);
        regionNlosses.back() += 1;
        break;
      }
      if (v[k].region != v[k - 1].region)
      {
        regionNlosses.back() += 1;
        // regionNlosses.emplace_back(k - prior);
        // prior = k;
        regionNlosses.emplace_back(0);
      }
      else if (v[k].year != v[k - 1].year)
      {
        regionNlosses.back() += 1;
      }
      k += 1;
    }
  }
  
  
  // return List::create(regionNlosses);
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  List aggLoss(Nregion), occLoss(Nregion);
  vec<std::pair<double*, double*>> rstPtr(Nregion);
  for (int i = 0; i < Nregion; ++i)
  {
    NumericVector x(regionNlosses[i]);
    rstPtr[i].first = &x[0];
    NumericVector y(regionNlosses[i]);
    rstPtr[i].second = &y[0];
    aggLoss[i] = x;
    occLoss[i] = y;
  }
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  if constexpr (true)
  {
    int k = 0;
    std::size_t j = 1, end = v.size(), t = 0; // t is pointer to the last element.
    double sumLoss = v[0].loss, maxLoss = v[0].loss;
    while (true)
    {
      if (j >= end)
      {
        rstPtr[k].first[t] = sumLoss;
        rstPtr[k].second[t] = maxLoss;
        // deqV[k].emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        // deqV.back().emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        break;
      }
      if (  v[j].region != v[j - 1].region )
      {
        // deqV[k].emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        // deqV.back().emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        rstPtr[k].first[t] = sumLoss;
        rstPtr[k].second[t] = maxLoss;
        maxLoss = v[j].loss;
        sumLoss = v[j].loss;
        k += 1;
        t = 0;
        // deqV.emplace_back(std::deque<std::pair<float, float>>());
      }
      else if ( v[j].year != v[j - 1].year )
      {
        // deqV[k].emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        // deqV.back().emplace_back(std::pair<float, float>(sumLoss, maxLoss));
        rstPtr[k].first[t] = sumLoss;
        rstPtr[k].second[t] = maxLoss;
        t += 1;
        maxLoss = v[j].loss;
        sumLoss = v[j].loss;
      }
      else
      {
        maxLoss = std::max(maxLoss, v[j].loss);
        sumLoss += v[j].loss;
      }
      j += 1;
    }
    
    // Rcout << "Nregion = " << Nregion << "\n";
    // Rcout << "k = " << k << "\n";
  }
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  cp.parFor(0, Nregion, [&](std::size_t i, std::size_t t)->bool
  {
    auto cmp = [](const double &x, const double &y)->bool { return x > y; };
    std::sort(rstPtr[i].first, rstPtr[i].first + regionNlosses[i], cmp);
    std::sort(rstPtr[i].second, rstPtr[i].second + regionNlosses[i], cmp);
    return false;
  });
  
  
  // Rcout << timer.tok() << "\n";
  // timer.tik();
  
  
  return List::create(Named("aggLoss") = aggLoss, Named("occLoss") = occLoss);
}




// [[Rcpp::export]]
List getSubareaPxIndicesAndImageDim(
    NumericVector subAreaLon, NumericVector subAreaLat, double pxLen)
{
  double minLon = 1e300, minLat = 1e300, maxLon = -1e300, maxLat = -1e300;
  for (int i = 0, iend = subAreaLon.size(); i < iend; ++i)
  {
    minLon = std::min(subAreaLon[i], minLon);
    minLat = std::min(subAreaLat[i], minLat);
    maxLon = std::max(subAreaLon[i], maxLon);
    maxLat = std::max(subAreaLat[i], maxLat);
  }
  minLon -= pxLen * 0.2;
  minLat -= pxLen * 0.2;
  int nx = int((maxLon - minLon) / pxLen) + 1; // ncol
  int ny = int((maxLat - minLat) / pxLen) + 1; // nrow
  IntegerVector matInd(subAreaLon.size());
  for (int i = 0, iend = subAreaLon.size(); i < iend; ++i)
  {
    int Xind = int((subAreaLon[i] - minLon) / pxLen);
    int Yind = int((subAreaLat[i] - minLat) / pxLen);
    int whichCol = Xind, whichRow = ny - Yind - 1;
    matInd[i] = whichCol * ny + whichRow;
  }
  return List::create(Named("nrow") = ny, Named("ncol") = nx, 
                      Named("subAreaPixelIndices") = matInd);
}


// [[Rcpp::export]]
NumericMatrix makeImage(List eventLoss, List subareaPxIndicesAndImageDim)
{
  IntegerVector ind = eventLoss[0];
  NumericVector loss = eventLoss[1];
  int nrow = subareaPxIndicesAndImageDim["nrow"];
  int ncol = subareaPxIndicesAndImageDim["ncol"];
  IntegerVector subAreaPixelIndices = 
    subareaPxIndicesAndImageDim["subAreaPixelIndices"];
  NumericMatrix img(nrow, ncol);
  for (int i = 0, iend = ind.size(); i < iend; ++i)
  {
    img[subAreaPixelIndices[ind[i]]] += loss[i];
  }
    
  return img;
}


// furtherestSigma: how many standard deviation is center to origin.
// [[Rcpp::export]]
NumericMatrix makeGauKernel(int windowSize, double furtherestSigma = 2)
{
  int halfSize = windowSize / 2;
  int size = halfSize * 2 + 1;
  NumericMatrix rst(halfSize * 2 + 1, halfSize * 2 + 1);
  double cx = halfSize, cy = halfSize;
  double S = 0;
  double sd = std::sqrt(cx * cx + cy * cy) / 2;
  for (int i = 0; i < size; ++i)
  {
    for (int j = 0; j < size; ++j)
    {
      double d = ((i - cx) * (i - cx) + (j - cy) * (j - cy)) / (2 * sd * sd);
      rst(j, i) = std::exp(-d);
      S += rst(j, i);
    }
  }
  S = 1.0 / S;
  for (auto &x: rst) x *= S;
  return rst;
}



// Y have been initialized with 0.
void sparseConv(double *X, int nrow, int ncol, 
                int *nonzeroInd, int nonzeroIndSize,
                double *Y, double *kern, int kernSize, 
                bool zeroYfirst,
                Charlie::VecPool &vp)
{
  if (zeroYfirst) std::fill(Y, Y + nrow * ncol, 0.0);
  auto visited = vp.lend<char>(nrow * ncol);
  std::fill(visited.begin(), visited.end(), false);
  int halfKern = kernSize / 2; // kernSize is always odd.
  auto getRowCol = [nrow, ncol](int k, int &row, int &col)->void
  {
    col = k / nrow;
    row = k % nrow;
  };
  
  
  // Convolve nearby pixels of (row, col).
  auto conv = [X, ncol, nrow, kern, halfKern, kernSize](int u, int v)->double // u: col, v: row
  {
    double S = 0;
    int virtualColStart = u - halfKern, virtualRowStart = v - halfKern;
    for (int i = std::max(0, u - halfKern),
         iend = std::min(ncol, u + halfKern + 1); i < iend; ++i)
    {
      for (int j = std::max(0, v - halfKern),
           jend = std::min(nrow, v + halfKern + 1); j < jend; ++j)
      {
        S += X[i * nrow + j] * kern[(
          i - virtualColStart) * kernSize + (j - virtualRowStart)];
      }
    }
    return S;
  };
  
  
  for (int i = 0; i < nonzeroIndSize; ++i)
  {
    int k = nonzeroInd[i];
    int row, col; getRowCol(k, row, col);
    for (int u = std::max(0, col - halfKern), 
         uend = std::min(ncol, col + halfKern + 1); 
         u < uend; ++u)
    {
      for (int v = std::max(0, row - halfKern),
           vend = std::min(nrow, row + halfKern + 1);
           v < vend; ++v)
      {
        if (visited[u * nrow + v]) continue;
        visited[u * nrow + v] = true;
        Y[u * nrow + v] = conv(u, v);
      }
    }
  }
  
  
  vp.recall(visited);
}


// [[Rcpp::export]]
NumericMatrix testSparseConv(NumericMatrix X, NumericMatrix kern)
{
  NumericMatrix Y(X.nrow(), X.ncol());
  Charlie::VecPool vp;
  auto nonzeroInd = vp.lend<int> (X.size());
  nonzeroInd.resize(0);
  for (int i = 0, iend = X.size(); i < iend; ++i)
  {
    if (X[i] == 0) nonzeroInd.emplace_back(i);
  }
  sparseConv(&X[0], X.nrow(), X.ncol(), 
             &nonzeroInd[0], nonzeroInd.size(),
             &Y[0], &kern[0], kern.nrow(), 
             false, vp);
  vp.recall(nonzeroInd);
  return Y;
}




// [[Rcpp::export]]
NumericMatrix imgConv(NumericMatrix X, NumericMatrix kern)
{
  int kernSize = kern.nrow(), halfKern = kernSize / 2;
  int nrow = X.nrow(), ncol = X.ncol();
  double *K = &kern[0];
  auto conv = [&](int Xc, int Yc)->double
  {
    int colStart = Xc - halfKern, colEnd = Xc + halfKern + 1;
    int rowStart = Yc - halfKern, rowEnd = Yc + halfKern + 1;
    double S = 0;
    for (int i = std::max(0, colStart), iend = std::min(colEnd, ncol); 
         i < iend; ++i)
    {
      for (int j = std::max(0, rowStart), jend = std::min(rowEnd, nrow);
           j < jend; ++j)
      {
        S += K[(i - colStart) * kernSize + (j - rowStart)] * X[i * nrow + j];
      }
    }
    return S;
  };
  NumericMatrix Y(nrow, ncol);
  for (int i = 0; i < ncol; ++i)
  {
    for (int j = 0; j < nrow; ++j)
      Y[i * nrow + j] = conv(i, j);
  }
  return Y;
} 





/*
// function(poolEventLosses, lossesOfEventsToBeReplaced, 
//          K = 1000L, relativeDiff = TRUE)
// [[Rcpp::export]]
List  findCandidates001(NumericVector poolEventLosses, 
                    NumericVector lossesOfEventsToBeReplaced,
                    int K = 1000, bool relativeDiff = true, 
                    int maxCore = 10000)
{
  if (K > poolEventLosses.size()) stop("K is too high.");
  List rst(lossesOfEventsToBeReplaced.size());
  vec<int*> rstPtr(rst.size());
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    IntegerVector v(K);
    rst[i] = v;
    rstPtr[i] = &v[0];
  }
  Charlie::ThreadPool cp(std::move(maxCore));
  vec<vec<std::pair<double, int>> > cntrV(maxCore);
  const bool useRelativeErr = relativeDiff;
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    auto &cntr = cntrV[t];
    cntr.resize(poolEventLosses.size());
    double thisLoss = lossesOfEventsToBeReplaced[i];
    for (int k = 0, kend = poolEventLosses.size(); k < kend; ++k)
    {
      if (useRelativeErr) cntr[k].first = std::abs((
        poolEventLosses[k] - thisLoss) / (poolEventLosses[k] + thisLoss));
      else cntr[k].first = std::abs(poolEventLosses[k] - thisLoss);
      cntr[k].second = k;
    }
    std::partial_sort(cntr.begin(), cntr.begin() + K, cntr.end());
    for (int k = 0; k < K; ++k) rstPtr[i][k] = cntr[k].second;
    return false;
  };
  cp.parFor(0, lossesOfEventsToBeReplaced.size(), f);
  return rst;
}
*/









#undef vec



