# pragma once


#include "charlieThreadPool.hpp"
#include "miniPCG.hpp"


#ifndef vec
#define vec std::vector
#endif


#define RNG MiniPcg32
#include <random>

namespace Charlie {


// =============================================================================
// Learning rate schedule: LinearLRschedule or ExpLRschedule
// Does not apply to every optimization problem.
// =============================================================================
template<typename ing, typename num>
struct LinearLRschedule
{
  num current, end, delta;
  void set(num start, num end, ing NstepToEnd)
  {
    this->end = end;
    delta = (end - start) / NstepToEnd;
    current = start;
  }
  LinearLRschedule(num start, num end, ing NstepToEnd)
  {
    set(start, end, NstepToEnd);
  }
  LinearLRschedule() {}
  void update()
  {
    current = std::max(end, current + delta);
  }
  num generate() { return current; }  // Read only, thread safe.
};


template<typename ing, typename num>
struct ExpLRschedule
{
  num current, end, delta;
  ExpLRschedule(num start, num end, ing NstepToEnd):
    end(end)
  {
    delta = std::pow(end / start, 1.0 / NstepToEnd);
    current = start;
  }
  ExpLRschedule() {}
  num update() { current = std::max(end, current * delta); }
  num generate() { return current; } // Read only, thread safe.
};
// =============================================================================


template<typename ing, typename GA>
struct cmp
{
  GA *x;
  cmp() {}
  cmp(GA *x): x(x) {}
  bool operator() (ing i, ing j) { return x->getval(i) < x->getval(j); }
};


template<typename ing, typename num, typename GA>
inline void ParaReproduce(
    GA *G, RNG *rngs, ing *odr, ing *parent, Charlie::ThreadPool *tp)
{
  tp->parFor(0, G->popuSize() - G->survivalSize(), [&](
      std::size_t objI, std::size_t t)->bool
  {
    G->reproduceTo(parent[objI], odr[G->survivalSize() + objI], rngs[objI]);
    return false;    
  }, 1);
}




// ======================================================================================
// The GA class should have the following member functions:
// 1.  void run(ing i): compute the function value of the i_th candidate.

// 2.  num & getval(ing i): return reference to the current function value of the i_th 
//     candidate.

// 3.  void copyTo(ing i, ing j): copy the i_th candidate's content to the j_th candidate.

// 4.  void reproduceTo(ing i, ing j, RNG &rng): i_th candidate reproduce a child, 
//     placing it in the
//     container occupied by j. This amounts to (i) copyTo(i, j), (ii) add noises to 
//     candidate j's parameters, (iii) run(j).

// 5.  ing popuSize(): return population size.

// 6.  ing survivalSize(): return survival population size in each generation.

// 7.  void actionsBetweenGenerations(ing iter): chance to adjust global parameters in GA 
//     after the iter_th generation.

// 8.  vec<ing> * corder(): return a pointer to the container for storing 
//     indices of candidates ordered by their function values. The container in the
//     GA class does not need to be initialized.
// ======================================================================================
// ======================================================================================
// Minimization.
// reproduceSelection: deterministic or randomized selection of parents for reproduction.
// 0: deterministic scan. Let every parent produces one child sequentially, and cycle.
// 1: linear prob kernel. 
// 2: inverse kernel.
// Return a numeric vector of vectors. result[i][j] is the function value of the j_th
// best candidate in the i_th generation.
// ======================================================================================
template <typename ing, typename num, typename GA>
vec<vec<num> > runGAobj(
    GA &initializedG, std::string &&reproduceSelectionMethod,
    ing NcandidateToSaveLearningCurve,
    ing maxIter, ing randomSeed, 
    Charlie::ThreadPool *cp, bool verbose = true)
{
  ing reproduceSelection = 0;
  if (reproduceSelectionMethod == "linearProb") reproduceSelection = 1;
  else if (reproduceSelectionMethod == "invProb") reproduceSelection = 2;
  
  
  GA &G = initializedG;
  ing popuSize = G.popuSize();
  ing survivalSize = G.survivalSize();
  NcandidateToSaveLearningCurve = std::min(
    survivalSize, NcandidateToSaveLearningCurve);
  
  
  vec<ing> &objfodr = *G.corder();
  objfodr.resize(popuSize);
  std::iota(objfodr.begin(), objfodr.end(), 0);
  
  
  vec<num> pvec;
  if (reproduceSelection != 0)
  {
    pvec.resize(survivalSize);
    if (reproduceSelection == 1)
    {
      
      for (ing i = 0; i < survivalSize; ++i) 
        pvec[i] = survivalSize - i;
    }
    else
    {
      for (ing i = 0; i < survivalSize; ++i) 
        pvec[i] = 1.0 / (i + 1);
    }
    std::partial_sum(pvec.begin(), pvec.end(), pvec.begin());
    for (auto & x: pvec) x /= pvec.back();
    pvec.back() = 1.0;
  }
  
  
  vec<ing> parent(popuSize - survivalSize);
  for (ing i = 0, iend = parent.size(); i < iend; ++i) 
    parent[i] = i % survivalSize;
  
  
  auto U = std::uniform_real_distribution<num> (0, 1);
  vec<RNG> rngs(parent.size());
  for (ing i = 0, iend = rngs.size(); i < iend; ++i) 
    rngs[i].seed(randomSeed + 1 + i);
  
  
  vec<vec<num> > learningCurve(maxIter, vec<num> (NcandidateToSaveLearningCurve));
  
  
  std::partial_sort(objfodr.begin(), objfodr.begin() + G.survivalSize(), 
                    objfodr.end(), cmp<ing, GA>(&G));
  
  
  for (ing iter = 0; iter < maxIter; ++iter)
  { 
    if (verbose)
    {
#ifdef Rcpp_hpp 
      Rcpp::Rcout << "Lowest object function value = " << G.getval(objfodr[0]) << "\n";
#else
      std::cout << "Lowest object function value = " << G.getval(objfodr[0]) << "\n";
#endif
    }
    
    
    if (reproduceSelection != 0)
    {
      for (ing k = survivalSize; k < popuSize; ++k)
      {
        ing u = std::lower_bound(pvec.begin(), pvec.end(), U(rngs[0])) - pvec.begin();
        parent[k - survivalSize] = u;
      }
    }
    
    
    ParaReproduce<ing, num, GA> (
        &G, &rngs[0], &objfodr[0], &parent[0], cp);
    
    
    std::partial_sort(objfodr.begin(), objfodr.begin() + G.survivalSize(), 
                      objfodr.end(), cmp<ing, GA>(&G));
    
    
    for (ing k = 0; k < NcandidateToSaveLearningCurve; ++k)
      learningCurve[iter][k] = G.getval(objfodr[k]);
    
    
    G.actionsBetweenGenerations(iter);
  }
  
  
  std::sort(objfodr.begin(), objfodr.end(), cmp<ing, GA>(&G));
  
  if (verbose) 
  {
#ifdef Rcpp_hpp
    Rcpp::Rcout << "Lowest object function value = " << 
      G.getval(objfodr[0]) << "\n";
#else
    std::cout << "Lowest object function value = " << 
      G.getval(objfodr[0]) << "\n";
#endif
  }
  
  return learningCurve;
}




}








