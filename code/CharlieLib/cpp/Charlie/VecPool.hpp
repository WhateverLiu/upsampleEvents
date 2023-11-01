#pragma once
#include "isVector.hpp"
#include "mmcp.hpp"


namespace Charlie {


// =============================================================================
// VecPool is not thread-safe. Create a VecPool for each thread.
// You have done enough. There has been so many ideas. The most advanced 
//   attempt we have tried is (i) write a class vec that inherits from std::vector,
//   (ii) inside the constructor and destructor of this new class, we would do
//   std::this_thread::get_id() to find which VecPool should we use.
//   This way, we can use vec freely and exactly like how we
//   would use std::vector. However, std::this_thread::get_id() is not an integer
//   and is not an index. We would need to use a hashtable to find the VecPool
//   lendn the thread id. Attempt of inheriting std::thread will get hairy fast.
//   
//   Now, don't think too much about persistent vectors or temp vectors inside
//     a class. It would be only a matter of putting the VecPool reference in
//     the function argument, or better, putting it as a function member and 
//     writing a destructor if you want to be more diligent.
// =============================================================================
class VecPool
{
private:
  std::vector<std::vector<char>> X;
  
  
  // Will zero both vectors' sizes. Capacities do not change.
  template <typename S, typename T>
  void darkSwap(std::vector<S> &s, std::vector<T> &t)
  {
    static_assert(sizeof(s) == sizeof(t));
    static_assert(sizeof(char*) == sizeof(std::size_t));
    static_assert(sizeof(s) == 3 * sizeof(char*));
    s.resize(0);
    t.resize(0);
    if constexpr ( std::is_same<S, T>::value ) { s.swap(t); }
    else
    {
      if constexpr ( std::is_same<T, char>::value ) darkSwap(t, s);
      else
      {
        if constexpr (std::is_same<S, char>::value) // s is char vector, but t is not.
        {
          std::size_t memT[3];
          mmcp(memT, &s, sizeof(s));
          memT[2] -= (memT[2] - memT[0]) % sizeof(T);
          mmcp(&s,  &t,  sizeof(s));
          mmcp(&t, memT, sizeof(t));
        }
        else
        {
          std::size_t memT[3]; mmcp(memT, &s, sizeof(s));
          std::size_t memS[3]; mmcp(memS, &t, sizeof(t));
          memT[2] -= (memT[2] - memT[0]) % sizeof(T);
          memS[2] -= (memS[2] - memS[0]) % sizeof(S);
          mmcp(&s, memS, sizeof(s));
          mmcp(&t, memT, sizeof(t));
        }
      }
    }
  }
  
  
  /*
  // ===========================================================================
  // Swap a std::vector<T> and a std::vector<char>.
  // ===========================================================================
  template<typename T> // rst is of size 0.
  void T2char(std::vector<T> &x, std::vector<char> &rst)
  { 
    static_assert(sizeof(x) == sizeof(std::vector<char>));
    x.resize(0); // Call destructors on all elements but maintain capacity.
    char mem[sizeof(x)]; // Start swapping.
    mmcp(mem,  &rst, sizeof(x));
    mmcp(&rst, &x,   sizeof(x));
    mmcp(&x,   mem,  sizeof(x));
  } 
  
  
  // ===========================================================================
  // Swap a std::vector<char> and a std::vector<T>.
  // ===========================================================================
  template<typename T>
  std::vector<T> char2T(std::vector<char> &x)
  {  
    static_assert(sizeof(x) == sizeof(std::vector<T>));
    x.resize(0);
    std::vector<T> rst;
    std::size_t mem[3];
    mmcp(mem, &x, sizeof(x));
    // =========================================================================
    // Adjust the capacity pointer to ensure the capacity is a multiple 
    //   of sizeof(T).
    // =========================================================================
    mem[2] -= (mem[2] - mem[0]) % sizeof(T);
    mmcp(&x,  &rst, sizeof(x));
    mmcp(&rst, mem, sizeof(x));
    return rst;
  }
  */
  
  
  // ===========================================================================
  // Dispatch and recall containers.
  // ===========================================================================
  template<typename T>
  std::vector<T> lendCore(std::size_t size)
  {
    if (X.size() == 0)
    {
      std::vector<T> rst(size + 1); // Always prefer slightly larger capacity.
      rst.pop_back();
      return rst;
    }
    std::vector<T> rst;
    darkSwap(rst, X.back());
    // auto rst = char2T<T> (X.back());
    rst.resize(size + 1); // Always prefer slightly larger capacity.
    rst.pop_back();
    X.pop_back();
    return rst;
  }
  
  
  // ===========================================================================
  // Recall the vector of an arbitrary type. The function will first clear
  // the vector --- call desctructor on all elements, and then fetch the
  // container back to pool.
  // ===========================================================================
  template<typename T>
  void recallCore(std::vector<T> &x) // invalidates everything related to x.
  { 
    if (x.capacity() == 0) return;
    x.resize(0); // Does not deallocate the container. This MUST go first
    // because T might not be trivial.
    // =========================================================================
    // Do not directly do X.emplace_back(x). This will have no effect since
    //   T2char(x).size() == 0. emplace_back ignores zero-size vectors. Tricky!
    // =========================================================================
    X.emplace_back(std::vector<char>());
    darkSwap(x, X.back());
    // T2char(x, X.back());
  }
  

public:
  std::vector<std::vector<char>> & getPool() { return X; }
  void test0(); // Test if std::vector is implemented as 3 pointers.
  void test1(); // Test if vectors fetched from or recalled to VecPool are handled correctly.
  void test2(); // Test if nested vectors can be handled correctly.
  // ===========================================================================
  // Treat typename... as a single keyword.
  // Ellipsis on the left of a type definer means the object is a parameter
  //   pack, on the right means unpacking.
  // ===========================================================================
  template <typename T, typename... Args>
  auto lend(std::size_t size, Args... restSizes) // restSizes is a pack.
  {
    // =========================================================================
    // sizeof...() is a single operator. It lends the number of parameters 
    //   in restSizes.
    // =========================================================================
    if constexpr (sizeof...(restSizes) == 0) return lendCore<T> (size);
    // =========================================================================
    // The other branch must be wrapped inside else {} to let the compiler
    //   know they are mutually exclusive.
    // =========================================================================
    else 
    { 
      // =======================================================================
      // lend<T>(restSizes...) will take the first parameter in the unpacked as
      //   the size argument, and the rest as a subpack.
      // =======================================================================
      auto rst = lendCore<decltype( lend<T>(restSizes...) )> (size);
      for (auto &x: rst) lend<T>(restSizes...).swap(x);
      return rst;
    }
  }
  
  
  // ===========================================================================
  // Recall the vector of an arbitrary type. If the type is also a vector or
  // nested vector of any depth, the function will recursively fetch all 
  // containers inside `x` back to pool.
  // ===========================================================================
  template <typename T> 
  void recall(std::vector<T> &x)
  {
    if constexpr (isVector<T>()())
    { 
      for (auto i = x.rbegin(); i != x.rend(); ++i) recall(*i);  
      recallCore(x);
    }
    else recallCore(x);
  }
  
  
  // Exchange the input vector for another vector.
  template <typename Out, typename In, typename... Args>
  auto swap ( std::vector<In> &x, Args... restSizes )
  {
    // If the input vector is not nested and output vector is also not nested,
    //   take a shortcut.
    static_assert(sizeof...(restSizes) != 0);
    if constexpr ( !isVector<In>()() and !isVector<Out>()() and 
                     sizeof...(restSizes) == 1 )
    {
      std::vector<Out> rst;
      darkSwap(rst, x);
      rst.resize(restSizes...);
      return rst;
    }
    else
    {
      recall(x);
      return lend(restSizes...);
    }
  }
  
  
  VecPool()
  {
    std::srand(std::time(nullptr));
    test0();
    test1();
    test2();
    std::vector<std::vector<char>>().swap(X);
  }
  
  
};


void VecPool::test0()
{
  static_assert( sizeof(std::size_t) == sizeof(void*) );
  
  
  // ===========================================================================
  // Test if std::vector is implemented in the way that can be exploited.
  // ===========================================================================
  unsigned fullSize = std::rand() % 31 + 17;
  unsigned subSize = fullSize - std::rand() % 13;
  
  
  typedef std::tuple<char, char, char> tupe; // Arbitrarily selected type.
  std::vector<tupe> a(fullSize); // Just for creating a vector of an arbitrary type.
  a.resize(subSize); // Downsize the vector.
  
  
  // ===========================================================================
  // To make VecPool work, vector header must contain exactly 3 pointers, 
  //   e.g. p1, p2, p3.
  // ===========================================================================
  static_assert( sizeof(a) == sizeof(std::size_t) * 3 );
  auto ptr = (std::size_t*)(&a); // Read but won't write.
  bool sizePointerCool     = ptr[0] + sizeof(tupe) * subSize  == ptr[1];
  bool capacityPointerCool = ptr[0] + sizeof(tupe) * fullSize == ptr[2];
  
  
  if (!sizePointerCool or !capacityPointerCool) throw std::runtime_error(
      "VecPool initialization: Layout of std::vector header is unsupported.\n");
}


// ===========================================================================
// Test 1 to run during initialization. Check memory addresses and if the 
//   container dispatched from VecPool can call push()/emplace_back().
// ===========================================================================
void VecPool::test1()
{
  std::size_t asize = std::rand() % 7 + 3;
  auto a = lend<float> (asize);
  a.resize(asize);
  void *aptr = a.data();
  
  
  int bsize = std::rand() % 5 + 3;
  auto b = lend<std::tuple<short, short, short> >(bsize);
  void *bptr = b.data();
  
  
  recall(b);
  recall(a);
  
  
  if ( (void*)(X.back().data()) != aptr ) throw std::runtime_error(
    "VecPool initialization: Last in pool is a different container.\n");
  
  
  if ( (void*)(X[X.size() - 2].data()) != bptr ) throw std::runtime_error(
    "VecPool initiailization: Second last in pool is a different container.\n");
  
  
  typedef std::tuple<char, char, char, char, char, char, char> tupe;
  int lastContainerByteCapa = X.back().capacity();
  int pushTimes = lastContainerByteCapa / sizeof(tupe) + 5;
  auto c = lend<tupe>(0);
  
  
  // =========================================================================
  // This will lend segmentation fault if c's byte capacity is not a 
  //   multiple of sizeof(tupe).
  // =========================================================================
  for (int i = 0; i < pushTimes; ++i) c.emplace_back(tupe());
  int dsize = c.capacity() - c.size() + 3;
  
  
  auto d = lend<double>(dsize);
  auto dptr = d.data() + 1;
  recall(d);
  if (std::size_t(X.back().data()) + sizeof(double) != std::size_t(dptr) ) 
    throw std::runtime_error(
        "VecPool initiailization: Last in pool is a different container -- 2.\n");
  
  
  // ===========================================================================
  // Test swap().
  // ===========================================================================
  auto e = lend<double> (rand() % 11 + 9);
  for (auto &x: e) x = rand() * 1e-8;
  double S = std::accumulate(e.begin(), e.end(), 0.0);
  auto edata = e.data();
  auto f = swap<short> ( e, e.size() / 3 * 2 );
  double Sf = std::accumulate(f.begin(), f.end(), 0.0);
  if ((char*)f.data() != (char*)edata) throw std::runtime_error(
    ".swap() is not working correctly. Some garbage value = " + 
      std::to_string(Sf + S) + ".\n");
}


template <typename T>
void FindAddrsOfAllVectorsInNestedVector(T &x, std::vector<std::size_t> &addrs)
{
  if constexpr (!Charlie::isVector<T>()()) return;
  else
  {
    static_assert(sizeof(std::size_t) == sizeof(x.data()));
    addrs.emplace_back(std::size_t(x.data()));
    for (auto &u: x) FindAddrsOfAllVectorsInNestedVector(u, addrs);
  }
}


// ===========================================================================
// Test 2 to run during initialization. Test if VecPool correctly deals with
//   nested vectors.
// ===========================================================================
void VecPool::test2()
{
  int a[7];
  for (int i = 0; i < 7; ++i) a[i] = rand() % 3 + 1;
  auto x = lend<float> ( // x is a 10-d vector.
    a[0], a[1], a[2], a[3], a[4], a[5], a[6]);
  std::vector<std::size_t> addresses;
  FindAddrsOfAllVectorsInNestedVector(x, addresses);
  recall(x);
  for (int i = X.size() - 1, j = 0; i >= 0; --i, ++j)
  { 
    auto addr = std::size_t(X[i].data());
    if (addresses[j] != addr) throw std::runtime_error(
          "VecPool initialization: Nested vector test failed.\n");
  }
} 


// =============================================================================
// Code pattern goes like this:
// =============================================================================
// VecPool vp;
// auto a = vp.lend<int> (3);
// auto b = vp.lend<double>(5, 7);
// auto c = vp.lend<std::tuple<short, short, short> >(2);
// 
// ......
// 
// vp.recall(c);
// vp.recall(b);
// vp.recall(a);
// =============================================================================
// It is strongly recommended that the first dispatched gets recalled last,
//   so the next time the code segment will always the same containers.
// =============================================================================


// ===========================================================================
// Code pattern for creating a vector of vectors of different sizes.
// ===========================================================================
// auto v = vp.lend< vec<vec<int>> > (3); // A vector of vector of vectors.
// for (int i = 0; i < 3; ++i)
// {
//   vp.lend<vec<int>>(i * 2 + 1).swap(v[i]);
//   for (int j = 0, jend = v[i].size(); j < jend; ++j)
//     vp.lend<int>(j + 5).swap(v[i][j]);
// }




}