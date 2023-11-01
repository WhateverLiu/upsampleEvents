// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie.hpp"


#define vec std::vector


// LRS is learning rate schedule.
template<typename ing, typename num, typename LearningRateFun>
struct BananaFun
{
  ing Nsurvival;
  LearningRateFun learningR;
  vec<std::pair<num, num> > param;
  vec<num> funval;
  vec<ing> odr;
  
  
  vec<ing> * corder() { return &odr; }
  
  
  void run(ing k)
  {
    num val1 = (1 - param[k].first);
    num val2 = param[k].second - param[k].first * param[k].first;
    funval[k] = val1 * val1 + val2 * val2 * 100;
  }
  
  
  num & getval(ing i)
  {
    return funval[i];
  }
  
  
  void copyTo(ing i, ing j)
  {
    param[j] = param[i];
    funval[j] = funval[i];
  }
  
  
  ing popuSize() { return param.size(); }
  ing survivalSize() { return Nsurvival; }
  
  
  void adjustParam() // After every generation. Single thread.
  {
    learningR.update();
  }
  
  
  void reproduceTo(ing i, ing j, Charlie::MiniPCG &rng)
  {
    copyTo(i, j);
    num learningRate = learningR.generate();
    std::uniform_real_distribution<num> U(-learningRate, learningRate);
    param[j].first += U(rng);
    param[j].second += U(rng);
    run(j);
  }
  
  
  BananaFun(){}
  BananaFun(std::pair<num, num> initParam, 
            ing popuSize, ing survivalSize,
            num iniNoise, num minNoise,
            ing maxGen, ing NgenerationsTillMinLearningRate):
    Nsurvival(survivalSize)
  {
    param.assign(popuSize, initParam);
    funval.resize(popuSize);
    learningR.set(iniNoise, minNoise, NgenerationsTillMinLearningRate);
    run(0);
    std::fill(funval.begin() + 1, funval.end(), funval.front());
  }
  
  
  void actionsBetweenGenerations(ing iter){}
};




// [[Rcpp::export]]
std::vector<std::vector<double> > testGA(
    Rcpp::NumericVector initxy, 
    double initNoise, double minNoise,
    int popuSize, int survivalSize, 
    int maxGen, int Ngen2minNoise,
    std::string reproduceSelection, int randomSeed,
    int maxCore)
{
  
  typedef BananaFun<int, double, Charlie::LinearLRschedule<int, double> > GAobj;
  GAobj ga(std::pair<double, double> (initxy[0], initxy[1]), 
           popuSize, survivalSize, initNoise, minNoise, maxGen, Ngen2minNoise);
  
  
  // template <typename ing, typename num, typename GA>
  // runGAobj(
  //   GA &initializedG, std::string reproduceSelectionMethod,
  //   ing NcandidateToSaveLearningCurve,
  //   ing maxIter, ing randomSeed, 
  //   CharlieThreadPool *cp, bool verbose = true)
  
  
  int NcandidateToSaveLearningCurve = 10;
  
  
  Charlie::ThreadPool cp(std::move(maxCore));
  Charlie::VecPool vp;
  auto rst = Charlie::runGAobj<int, double, GAobj> (
    ga, std::move(reproduceSelection), NcandidateToSaveLearningCurve, 
    maxGen, randomSeed, cp, vp, true);
  
  
  return rst;
} 


#undef vec



