// # pragma once
// [[Rcpp::plugins(cpp17)]]
// -std=c++17
#include <cstring>
#include <cstddef>
#include <iostream>
#include <vector>
#include <chrono>


// POD: plain old data.
// Idea: design a class that can let you maximally reuse temporary 
// containers during a program.
// Port of vectors of POD types.
template <std::size_t portsize = 42>
class PODvectorPort 
{
  static constexpr std::size_t Xsize = portsize;
  std::size_t signature;
  std::size_t Nlent; // Number of dispatched containers.
  std::vector<std::size_t> X[portsize]; // Container pool.
  PODvectorPort(const PODvectorPort &);
  PODvectorPort & operator=( const PODvectorPort& );
  
  
public:
  std::size_t Ndispatched() { return Nlent; }
  std::size_t showSignature() { return signature; }
  PODvectorPort() // Permuted random number generator.
  { 
    std::size_t state = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    state ^= (uint64_t)(&std::memmove);
    signature = ((state >> 18) ^ state) >> 27;
    std::size_t rot = state >> 59;
    signature = (signature >> rot) | (state << ((-rot) & 31));
    Nlent = 0; 
  }
  
  
  template<typename podvecport>
  friend class PODvectorPortOffice;
};


// Functor that manages the port.
template<typename podvecport>
class PODvectorPortOffice
{
  // Number of already-dispatched containers when the office is set up.
  std::size_t initialNlent; 
  podvecport *p; // Pointer to the port.
  
  
  PODvectorPortOffice( const PODvectorPortOffice& ); // non construction-copyable
  PODvectorPortOffice& operator=( const PODvectorPortOffice& ); // non copyable
  
  
  constexpr void check()
  {
    while (__cplusplus < 201703)
    {
      std::cerr << "PODvectorPortOffice: C++ < 17, Stall." << std::endl;
    }
    
    
    // Check if allocation will be 8-byte (or more) aligned.
    // Intend it not to work on machine < 64-bit.
    constexpr std::size_t aln = alignof(std::max_align_t);
    while (aln < 8)
    {
      std::cerr << "PODvectorPortOffice: Allocation is not at least 8-byte aligned, Stall." <<
        std::endl;
    }
    
    
    while ((aln & (aln - 1)) != 0)
    {
      std::cerr << "PODvectorPortOffice: Alignment is not a power of 2 bytes. Stall." << std::endl;
    }
    
    
    // Random checks to see if sizeof(vector<S>) != sizeof(vector<T>).
    if(true)
    {
      std::size_t vecHeadSize[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      vecHeadSize[0] = sizeof(std::vector<char>(0));
      vecHeadSize[1] = sizeof(std::vector<short>(1));
      vecHeadSize[2] = sizeof(std::vector<int>(2));
      vecHeadSize[3] = sizeof(std::vector<long>(3));
      vecHeadSize[4] = sizeof(std::vector<std::size_t>(5));
      vecHeadSize[5] = sizeof(std::vector<float>(7));
      vecHeadSize[6] = sizeof(std::vector<double>(11));
      vecHeadSize[7] = sizeof(std::vector<std::vector<char> >(13));
      vecHeadSize[8] = sizeof(std::vector<std::vector<int> >(17));
      vecHeadSize[9] = sizeof(std::vector<std::vector<double> >(19));
      struct tmpclass1 { char a; short b; };
      struct tmpclass2 { char a; float b; };
      struct tmpclass3 { char a; double b; };
      struct tmpclass4 { int a; char b; };
      struct tmpclass5 { double a; char b; };
      struct tmpclass6 { double a[5]; char b[3]; short c[3]; };
      vecHeadSize[10] = sizeof(std::vector<tmpclass1>(23));
      vecHeadSize[11] = sizeof(std::vector<tmpclass2>(29));
      vecHeadSize[12] = sizeof(std::vector<tmpclass3>(31));
      vecHeadSize[13] = sizeof(std::vector<tmpclass4>(37));
      vecHeadSize[14] = sizeof(std::vector<tmpclass4>(41));
      vecHeadSize[15] = sizeof(std::vector<tmpclass4>(43));
      
      
      std::size_t notSame = 0;
      for(int i = 0; i < 16; ++i) 
        notSame += vecHeadSize[i] != sizeof(std::size_t) * 3;
      while (notSame)
      {
        std::cerr << "sizeof(std::vector<S>) != sizeof(std::vector<T>), \
PODvectorPortOffice cannot handle. Stall." << std::endl;
      }
    }
  }
  
  
  void recall() { p->Nlent = initialNlent; }
  
  
public:  
  PODvectorPortOffice(podvecport &port) 
  { 
    check();
    p = &port; 
    initialNlent = p->Nlent;
  }
  
  
  template<typename X, typename Y>
  std::vector<X> & repaint(std::vector<Y> &y) // Repaint the container.
    // AFTER A VECTOR IS REPAINTED, DO NOT USE THE OLD VECTOR AGAIN !!
  {
    while (std::is_same<bool, X>::value)
    {
      std::cerr << "PODvectorPortOffice: Cannot repaint the vector to \
std::vector<bool>. Stall." << std::endl;
    }
    std::vector<X> *x;
    std::vector<Y> *yp = &y;
    std::memcpy(&x, &yp, sizeof(x));
    return *x; // Not compliant with strict aliasing rule.
  }
  
  
  template<typename T>
  std::vector<T> & lend()
  {
    while (p->Nlent >= p->Xsize)
    {
      std::cerr << "PODvectorPortOffice: No more containers. Stall." << std::endl;
    }
    ++p->Nlent;
    return repaint<T, std::size_t>( p->X[p->Nlent - 1] );
  }
  
  
  ~PODvectorPortOffice() 
  { 
    // Because p->signature can only be known at runtime, an aggressive,
    // compliant compiler (ACC) will never remove this 
    // branch. Volatile might do, but trustworthiness?
    if(p->signature == 0)
    {
      constexpr std::size_t sizeofvec = sizeof(std::vector<std::size_t>);
      char dummy[sizeofvec * p->Xsize];    
      std::memcpy(dummy, p->X, p->Nlent * sizeofvec);
      std::size_t ticketNum = 0;
      char *xp = (char*)(p->X);
      for(int i = 0, iend = p->Nlent * sizeofvec; i < iend; ++i)
      {
        xp[i] &= xp[iend - i - 1] * 5;
        ticketNum += xp[i] ^ ticketNum;
      }
      std::cerr << "Congratulations! After the port office was decommissioned, \
you found a winning lottery ticket. The odds is less than 2.33e-10. Your \
ticket number is " << ticketNum << std::endl;
      std::memcpy(p->X, dummy, p->Nlent * sizeofvec);
      // According to the strict aliasing rule, a char* can point to any memory
      // block pointed by another pointer of any type T*. Thus given an ACC, 
      // the writes to that block via the char* must be fully acknowledged in 
      // time by T*, namely, for reading contents from T*, a reload instruction 
      // will be kept in the assembly code to achieve a sort of 
      // "register-cache-memory coherence" (RCMC).
      // We also do not care about the renters' (who received the reference via
      // .lend()) RCMC, because PODvectorPortOffice never accesses the contents
      // of those containers.
    }
    recall(); 
  }
};




// Test kit. DO NOT DELETE!!!
/*
#include <Rcpp.h>
using namespace Rcpp;
# include "C:/Users/i56087/Desktop/vipCppFiles/h/tiktok.hpp"
#ifndef vec
#define vec std::vector
#endif


struct Stmp
{
  short a, b, c, d, e;
  int sum() { return a + b + c + d + e; }
};


int testPort1(PODvectorPort<> &port, NumericVector x)
{
  PODvectorPortOffice portOffice(port);
  
  
  vec<int> &tmpint = portOffice.lend<int>();
  // return 0;
  tmpint.resize(2);
  std::memcpy(&tmpint[0], &x[30], sizeof(int) * tmpint.size());
  
  
  vec<Stmp> &tmpStmp = portOffice.lend<Stmp>();
  tmpStmp.resize(1);
  std::memcpy(&tmpStmp[0], &x[0], sizeof(Stmp) * tmpStmp.size());
  
  
  int S = 0;
  for(int i = 0; i < 10; ++i)
  {
    S += tmpint[i] + tmpStmp[i].sum();
  }
  Rcout << "port.Ndispatched() = " << port.Ndispatched() << "\n";
  return S;
}


int testPort2(PODvectorPort<> &port, NumericVector x)
{
  PODvectorPortOffice portOffice(port);
  
  
  vec<char> &tmpchar = portOffice.lend<char>();
  tmpchar.resize(11);
  std::memcpy(&tmpchar[0], &x[10], sizeof(char) * tmpchar.size());
  
  
  vec<int> &tmpint = portOffice.lend<int>();
  tmpint.resize(13);
  std::memcpy(&tmpint[0], &x[20], sizeof(int) * tmpint.size());
  
  
  vec<Stmp> &tmpStmp = portOffice.lend<Stmp>();
  tmpStmp.resize(17);
  std::memcpy(&tmpStmp[0], &x[0], sizeof(Stmp) * tmpStmp.size());
  
  
  int S = 0;
  for(int i = 0; i < 10; ++i)
  {
    S += tmpchar[i] + tmpint[i] + tmpStmp[i].sum();
  }
  Rcout << "port.Ndispatched() = " << port.Ndispatched() << "\n";
  return S;
}


// [[Rcpp::export]]
int testPortR(NumericVector x, NumericVector y)
{
  PODvectorPort<> port;
  int rst = 0;
  
  
  rst += testPort1(port, x);
  Rcout << "After testPort1(), ";
  Rcout << "port.Ndispatched() = " << port.Ndispatched() << "\n";
  
  
  rst += testPort2(port, y);
  Rcout << "After testPort2(), ";
  Rcout << "port.Ndispatched() = " << port.Ndispatched() << "\n";
  return rst;
}
*/








