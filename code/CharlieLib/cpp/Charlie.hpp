// [[Rcpp::depends(RcppEigen)]]
#include "Charlie/LBFGSB.hpp"
#include "Charlie/tiktok.hpp"
#include "Charlie/mmcp.hpp"
#include "Charlie/isVector.hpp"
#include "Charlie/DarkSwapVec.hpp"
#include "Charlie/VecPool.hpp"
#include "Charlie/MiniPCG.hpp"
#include "Charlie/ThreadPool.hpp"
#include "Charlie/Sort.hpp"
#include "Charlie/LongestMonoSubseq.hpp"
#include "Charlie/MonoLinearInterpo.hpp"
#include "Charlie/GA.hpp"
#include "Charlie/xxhash.hpp"
// =============================================================================
// TempAlloc.hpp: do not attempt. No good! Too much trouble to use and to setup.
// Consider overload new and delete. But No! That will not allow you use
// std::vector any more! Also, the power of 2 expansion could easily
// deplete the memory! Stick to VecPool!
// =============================================================================
// #include "Charlie/TempAlloc.hpp"

