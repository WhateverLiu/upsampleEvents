#include "lbfgs-okazaki-yixuan/LBFGS.h"
#include "lbfgs-okazaki-yixuan/LBFGSB.h"


// =============================================================================
// Function object MUST have 1 method:
//   
//
// double operator()(double *x, double *&grad, int dim);
//    Return function value and assign grad to nullptr if the function does not
//    compute the gradient, otherwise it should write the gradient to space
//    pointed by grad.
//
//
// `hgrad` is the scaler for computing gradients via finite 
// differencing if gradient is not given. In that case, the step size h for 
// computing the gradient would be max(1, |x|) * hgrad.
// =============================================================================
namespace Charlie {


// L-BFGS and L-BFGS-B combo
struct LBFGSB 
{
  // Result parameters will be written into the initialization vector.
  // Return min function value and the number of iterations taken.
  template <typename objFun>
  std::pair<double, int> operator() (
      double *x,
      int dim, 
      objFun &f, 
      int max_iterations = 100,
      double *LB = nullptr, // Use L-BFGS if both LB and UB are nullptr, 
      double *UB = nullptr, //  otherwise use L-BFGS-B.
      double hgrad = 0,
      bool centralDiff = true,
      // =======================================================================
      // Read https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html
      // to understand parameter setting in LBFGSB. Do not change them unless
      // necessary.
      // =======================================================================
      double epsilon = 1e-5,
      double epsilon_rel = 1e-5,
      int m = 6, // History
      int past = 1,
      double delta = 1e-10,
      int max_submin     = 10,
      int max_linesearch = 20,
      double min_step = 1e-20,
      double max_step = 1e+20,
      double ftol = 1e-4,
      double wolfe = 0.9
  // =======================================================================
  )
  { 
    
    if (centralDiff) return run<objFun, true>(
        x,
        dim, 
        f, 
        LB,
        UB,
        hgrad,
        max_iterations,
        epsilon,
        epsilon_rel,
        m,
        past,
        delta,
        max_submin,
        max_linesearch,
        min_step,
        max_step,
        ftol,
        wolfe);
    
    
    return run<objFun, false>(
        x,
        dim, 
        f, 
        LB,
        UB,
        hgrad,
        max_iterations,
        epsilon,
        epsilon_rel,
        m,
        past,
        delta,
        max_submin,
        max_linesearch,
        min_step,
        max_step,
        ftol,
        wolfe);
  }
  
  
  template <typename objFun, bool centralDiff>
  std::pair<double, int> run (
      double *x, // Result parameters will be written into the initialization vector.
      int dim, 
      objFun &f, 
      double *LB = nullptr, 
      double *UB = nullptr,
      double hgrad = 0,
      // =======================================================================
      // Read https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html
      // to understand parameter setting in LBFGSB. Do not change them unless
      // necessary.
      // =======================================================================
      int max_iterations = 100,
      double epsilon = 1e-5,
      double epsilon_rel = 1e-5,
      int m = 6,
      int past = 1,
      double delta = 1e-10,
      int max_submin     = 10,
      int max_linesearch = 20,
      double min_step = 1e-20,
      double max_step = 1e+20,
      double ftol = 1e-4,
      double wolfe = 0.9
  // =======================================================================
  )
  {
    
    // Amount recommended in paper. MachinePrecision ^ (1/3) and 
    //   MachinePrecision ^ (1/2).
    if (centralDiff) hgrad = hgrad == 0 ? 6.05545445239334e-06 : hgrad;
    else hgrad = hgrad == 0 ? 1.49011611938477e-08 : hgrad;
    
    
    const bool unbounded = LB == nullptr and UB == nullptr;
    
    
    if (!unbounded) param.reset(
        max_iterations,
        m,
        epsilon,
        epsilon_rel,
        past,
        delta,
        max_submin,
        max_linesearch,
        min_step,
        max_step,
        ftol,
        wolfe
    ); else paramUnbounded.reset(
        max_iterations,
        m,
        epsilon,
        epsilon_rel,
        past,
        delta,
        max_linesearch,
        min_step,
        max_step,
        ftol,
        wolfe
    );
    
    
    if (!unbounded)
    {
      if (LB == nullptr)
      {
        lb.resize(dim);
        std::fill(&lb[0], &lb[0] + dim, -1e300);
      }
      else lb = Eigen::Map<Eigen::VectorXd> (LB, dim);
      
      
      if (UB == nullptr)
      {
        ub.resize(dim);
        std::fill(&ub[0], &ub[0] + dim, 1e300);
      }
      else ub = Eigen::Map<Eigen::VectorXd> (UB, dim);
    }
    
    
    // Stop thinking about micro-optimizing these copies! 
    // The speed and space cost is a tiny fraction of the internal 
    // library's if you read it!
    xv = Eigen::Map<Eigen::VectorXd> (x , dim);
    
    
    objFun *fptr = &f;
    auto fun = [dim, fptr, hgrad](
      Eigen::VectorXd &x, Eigen::VectorXd &grad)->double
      {
        double *g = &grad[0];
        double fval = (*fptr)(&x[0], g, dim);
        if (g != nullptr) return fval;
        for (int i = 0; i < dim; ++i) // Finite difference for gradient.
        { 
          double safe = x[i];
          double dx = std::abs(x[i]) * hgrad;
          x[i] += dx;
          double fhigh = (*fptr)(&x[0], g, dim);
          if (!centralDiff) grad[i] = (fhigh - fval) / dx;
          else
          { 
            x[i] = safe - dx;
            double flow = (*fptr)(&x[0], g, dim);
            grad[i] = (fhigh - flow) / (2 * dx);
          } 
          x[i] = safe;
        } 
        return fval;
      };
    
    
    double finalFval = 0;
    int totIter = 0;
    
    
    try {
      if (!unbounded)
      {
        totIter = solver.minimize(fun, xv, finalFval, lb, ub);
      }
      else
      {
        totIter = solverUnbounded.minimize(fun, xv, finalFval);
      }
    }
    catch( ... )
    {
      std::fill(&xv[0], &xv[0] + dim, 1e300);
      totIter = 0;
      finalFval = 1e300;
    }
    
    
    std::copy(&xv[0], &xv[0] + dim, x);
    return std::pair<double, int> (finalFval, totIter);
  }
  
  
  LBFGSpp::LBFGSBParam<double> param;
  LBFGSpp::LBFGSParam<double> paramUnbounded;
  LBFGSpp::LBFGSBSolver<double> solver;
  LBFGSpp::LBFGSSolver<double> solverUnbounded;
  Eigen::VectorXd lb, ub, xv;
  LBFGSB(): solver(param), solverUnbounded(paramUnbounded) {}
  
  
};



}




