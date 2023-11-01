#pragma once
// #include "LongestMonoSubseq.hpp"


namespace Charlie {
#define vec std::vector


// =============================================================================
// lendn y[size] and SORTED x[size], find the longest monotonic 
//   subsequence of y[] and the associated subsequence of x. Linearly 
//   interpolate a function y = f(x).
// The class should overload operator() for query.
//
// monoType == 0: increasing subsequence.
// monoType == 1: nondecreasing subsequence.
// monoType == 2: decreasing subsequence.
// monoType == 3: nonincreasing subsequence.
// =============================================================================
template <typename ing, typename numX, typename numY>
struct MonoLinearInterpo
{
  LongestMonoSubseq<ing, numY> lms;
  vec<numX> xsubseq;
  vec<numY> ysubseq;
  VecPool *vp;
  
  
  // Compute the longest subsequences of both x and y. Get ready for query().
  template<int monoType>
  void reset(numX *x, numX *xend, numY *y, VecPool &vp)
  {
    auto indx = vp.lend<ing>(xend - x);
    auto indxEnd = lms.template operator()<monoType> (
      y, y + (xend - x), indx.data(), vp);
    indx.resize(indxEnd - indx.data());
    
    
    if (xsubseq.capacity() == 0) vp.lend<numX>(indx.size()).swap(xsubseq);
    else xsubseq.resize(indx.size());
    if (ysubseq.capacity() == 0) vp.lend<numY>(indx.size()).swap(ysubseq);
    else ysubseq.resize(indx.size());
    
    
    for (int64_t i = 0, iend = indx.size(); i < iend; ++i)
    {
      xsubseq[i] = x[indx[i]];
      ysubseq[i] = y[indx[i]];
    }
    
    
    vp.recall(indx);
    this->vp = &vp;
  } 
  
  
  ~MonoLinearInterpo()
  {
    vp->recall(ysubseq);
    vp->recall(xsubseq);
  }
  
  
  numY operator()(numX &&x)
  {
    auto it = std::lower_bound(xsubseq.begin(), xsubseq.end(), x);
    if (it >= xsubseq.end()) return ysubseq.back();
    if (it <= xsubseq.begin()) return ysubseq.front();
    auto xlow = *(it - 1), xhigh = *it;
    double w = (x - xlow) / (xhigh - xlow);
    auto yi = it - xsubseq.begin();
    return (1 - w) * ysubseq[yi - 1] + w * ysubseq[yi];
  }
  
  
  numY operator()(numX &x) { return (*this)(std::move(x));  }
  
  
  void operator()(numX *x, numX *xend, numY *y)
  {
    for (int64_t i = 0, iend = xend - x; i < iend; ++i)
      y[i] = (*this)(x[i]);
  }
};


#undef vec
}




