#pragma once
#include "ThreadPool.hpp"
#include "VecPool.hpp"
// 
// 
// #ifndef vec 
// #define vec std::vector
// #endif


// Use Multithreaded merge sort if the vector is big enough and a pointer to
// an existing CharlieThreadPool is lendn.


namespace Charlie {


// Can only sort elements in a block of memory.
struct Sort
{

  template <typename num, typename Compare>  
  void operator() (
      num *x, num *xend, Compare &&f, 
      ThreadPool *cp = nullptr, // Using pointer not reference to imply if multithreading is on.
      VecPool *vp = nullptr)
  {
    if (cp == nullptr or xend - x < 3000 or cp->maxCore <= 1)
    { 
      std::sort(x, xend, f); return; 
    }
    
    
    uint64_t Nblock = std::round(std::exp2(uint64_t(
      std::round(std::log2(cp->maxCore * 4.0) / 2)) * 2));
    
    
    std::vector<uint64_t> offset;
    if (vp == nullptr) offset.resize(Nblock + 1);
    else vp->lend<uint64_t>(Nblock + 1).swap(offset);
    
    
    offset[0] = 0;
    double approxBlockSize = (xend - x + 0.0) / Nblock;
    for (int64_t i = 1, iend = offset.size() - 1; i < iend; ++i)
      offset[i] = uint64_t(std::round(approxBlockSize * i));
    offset.back() = xend - x;
    cp->parFor(0, Nblock, [&](std::size_t i, std::size_t t)->bool
    {
      std::sort(x + offset[i], x + offset[i + 1], f);
      return false;
    }, 1);
    
    
    std::vector<num> v;
    if (vp != nullptr) vp->lend<num>(xend - x).swap(v);
    else v.resize(xend - x);
    
    
    num *y = v.data();
    for (uint64_t delta = 1; delta < Nblock; delta *= 2)
    {
      cp->parFor(0, Nblock / (delta * 2), [&](std::size_t i, std::size_t t)->bool
      {
        std::size_t i2delta = i * 2 * delta;
        uint64_t &begin1 = offset[i2delta];
        uint64_t &end1   = offset[i2delta + delta];
        uint64_t &begin2 = end1;
        uint64_t &end2   = offset[(i + 1) * 2 * delta];
        std::merge(x + begin1, x + end1, x + begin2, x + end2, y + begin1, f);
        return false;
      }, 1);
      std::swap(x, y);
    }
    
    
    if (vp != nullptr) { vp->recall(v); vp->recall(offset); }
  }
  
  
  template <typename num>
  void operator() (num *x, num *xend,
                Charlie::ThreadPool *cp = nullptr,
                Charlie::VecPool *vp = nullptr)
  {
    auto cmp = [](const num &u, const num &v)->bool { return u < v; };
    (*this)(x, xend, cmp, cp);
  }
  
  
  template <typename Iter, typename Compare>  
  void operator() (Iter x, Iter xend, Compare &&f, 
                Charlie::ThreadPool *cp = nullptr,
                Charlie::VecPool *vp = nullptr)
  {
    (*this)(&*x, &*xend, f, cp);
  }
  
  
  template <typename Iter>
  void operator() (Iter x, Iter xend, 
                Charlie::ThreadPool *cp = nullptr,
                Charlie::VecPool *vp = nullptr)
  {
    (*this)(&*x, &*xend, cp);
  }
};


}










