// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../Charlie.hpp"


// [[Rcpp::export]]
void test1(int asize = 3) 
{
  Charlie::VecPool vp;
  auto a = vp.lend<float>(asize);
  a.resize(asize);
  void *aptr = a.data();
  
  
  auto b = vp.lend<std::tuple<short, short, short> >(asize + 2);
  void *bptr = b.data();
  
  
  vp.recall(b);
  vp.recall(a);
   
  
  if ( (void*)(vp.getPool().back().data()) != aptr )
    throw std::runtime_error("Last in pool is a different container.");
  
  
  if ( (void*)(vp.getPool()[vp.getPool().size() - 2].data()) != bptr )
    throw std::runtime_error("Second last in pool is a different container.");
  
  
  typedef std::tuple<char, char, char, char, char, char, char> tupe;
  int lastContainerByteCapa = vp.getPool().back().capacity();
  Rcout << "lastContainerByteCapa = " << lastContainerByteCapa << "\n";
  int pushTimes = lastContainerByteCapa / sizeof(tupe) + 5;
  auto c = vp.lend<tupe>(0);
  Rcout << "c.size() = " << c.size() << "\n";
  Rcout << "c.capacity() = " << c.capacity() << "\n";
  
  
  // This will lend segmentation fault if c's byte capacity is not a 
  //   multiple of sizeof(tupe).
  for (int i = 0; i < pushTimes; ++i)
    c.emplace_back(tupe());
  
  
  Rcout << "c.size() = " << c.size() << "\n";
  Rcout << "c.capacity() = " << c.capacity() << "\n";
  Rcout << "Capacity works correctly.\n";
  
  
  auto *cptr = c.data();
  vp.recall(c);
  if ( (void*)(vp.getPool().back().data()) != cptr )
    throw std::runtime_error("Last in pool is a different container -- 2.");
} 




struct TmpFun
{
  std::vector<std::size_t> *addrsPtr;
  TmpFun(std::vector<std::size_t> &addrs) { addrsPtr = &addrs; }
  
  template <typename T>
  void operator()(T &x) // Function will find all container's addresses.
  {
    if constexpr (!Charlie::isVector<T>()()) return;
    else
    {
      addrsPtr->emplace_back(std::size_t(x.data()));
      for (auto &u: x) (*this)(u);
    }
  }
};




// [[Rcpp::export]]
void test2(int seed)
{
  
  Charlie::MiniPCG rng(seed);
  Charlie::VecPool vp;
  auto U = std::uniform_int_distribution<int> (1, 5); // At most 5^10 elements could be allocated.
  int a[10];
  for (int i = 0; i < 10; ++i) a[i] = U(rng);
  auto x = vp.lend<float> ( // x is a 10-d vector.
    a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]);
  std::vector<std::size_t> addresses;
  auto f = TmpFun(addresses); f(x);
  
  
  vp.recall(x);
  std::size_t mib = 0;
  Rcout << "Total number of vectors allocated = " << vp.getPool().size() << "\n";
  for (int i = vp.getPool().size() - 1, j = 0; i >= 0; --i, ++j)
  {
    mib += vp.getPool()[i].capacity();
    auto addr = std::size_t(vp.getPool()[i].data());
    if (addresses[j] != addr)
    {
      Rcout << "Address when it's initialized = " << addresses[j] << ", ";
      Rcout << "Address back in pool = " << addr << "\n";
    }
  }
  Rcout << "Total number of MiB allocated = " << mib / (1024.0 * 1024) << "\n";
}






