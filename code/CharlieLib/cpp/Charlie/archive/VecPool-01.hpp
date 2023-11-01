namespace Charlie {


// =============================================================================
// VecPool is not thread-safe. Create a VecPool for each thread.
// =============================================================================
template <typename T>
struct mvec: public std::vector<T>
{
  thread_local static VecPool *pool;
  void make()
  {
    
  }
  mvec(std::size_t size) { resize(size); }
  mvec(std::size_t size, T &&t) { assign }
  
  ~mvec()
  {
    pool
  }
};




class VecPool
{
private:
  std::vector<std::vector<char>> X;
  
  
  // ===========================================================================
  // Swap a std::vector<T> and a std::vector<char>.
  // ===========================================================================
  template<typename T> // rst is of size 0.
  // void T2char(std::vector<T> &x, std::vector<char> &rst)
  void T2char(mvec<T> &x, std::vector<char> &rst)
  { 
    static_assert(sizeof(x) == sizeof(rst));
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
  // std::vector<T> char2T(std::vector<char> &x)
  void char2T(std::vector<char> &x, mvec<T> &rst)
  {  
    static_assert(sizeof(x) == sizeof(rst));
    x.resize(0);
    rst.resize(0);
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
  
  
  // ===========================================================================
  // Dispatch and recall containers.
  // ===========================================================================
  template<typename T>
  // std::vector<T> giveCore(std::size_t size)
  mvec<T> giveCore(std::size_t size)
  {
    if (X.size() == 0) 
    {
      mvec<T> rst(size + 1); // Always prefer slightly larger capacity.
      rst.pop_back();
      return rst;
    }
    auto rst = char2T<T> (X.back());
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
    x.resize(0); // Does not deallocate the container.
    // =========================================================================
    // Do not directly do X.emplace_back(rst). This will have no effect since
    //   T2char(x).size() == 0. emplace_back ignores zero-size vectors. Tricky!
    // =========================================================================
    X.emplace_back(std::vector<char>(0));
    T2char(x, X.back());
  }
  

public:
  std::vector<std::vector<char>> &getPool() { return X; }
  void test0(); // Test if std::vector is implemented as 3 pointers.
  void test1(); // Test if vectors fetched from or recalled to VecPool are handled correctly.
  void test2(); // Test if nested vectors can be handled correctly.
  // ===========================================================================
  // Treat typename... as a single keyword.
  // Ellipsis on the left of a type definer means the object is a parameter
  //   pack, on the right means unpacking.
  // ===========================================================================
  template <typename T, typename... Args>
  auto give(std::size_t size, Args... restSizes) // restSizes is a pack.
  {
    // =========================================================================
    // sizeof...() is a single operator. It gives the number of parameters 
    //   in restSizes.
    // =========================================================================
    if constexpr (sizeof...(restSizes) == 0) return giveCore<T> (size);
    // =========================================================================
    // The other branch must be wrapped inside else {} to let the compiler
    //   know they are mutually exclusive.
    // =========================================================================
    else 
    { 
      // =======================================================================
      // give<T>(restSizes...) will take the first parameter in the unpacked as
      //   the size argument, and the rest as a subpack.
      // =======================================================================
      auto rst = giveCore<decltype( give<T>(restSizes...) )> (size);
      for (auto &x: rst) give<T>(restSizes...).swap(x);
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
  
  
  VecPool()
  {
    std::srand(std::time(nullptr));
    test0();
    test1();
    test2();
    std::vector<std::vector<char>>(0).swap(X);
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
  auto a = give<float> (asize);
  a.resize(asize);
  void *aptr = a.data();
  
  
  int bsize = std::rand() % 5 + 3;
  auto b = give<std::tuple<short, short, short> >(bsize);
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
  auto c = give<tupe>(0);
  
  
  // =========================================================================
  // This will give segmentation fault if c's byte capacity is not a 
  //   multiple of sizeof(tupe).
  // =========================================================================
  for (int i = 0; i < pushTimes; ++i) c.emplace_back(tupe());
  int dsize = c.capacity() - c.size() + 3;
  
  
  auto d = give<double>(dsize);
  auto dptr = d.data() + 1;
  recall(d);
  if (std::size_t(X.back().data()) + sizeof(double) != std::size_t(dptr) ) 
    throw std::runtime_error(
        "VecPool initiailization: Last in pool is a different container -- 2.\n");
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
  auto x = give<float> ( // x is a 10-d vector.
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
// auto a = vp.give<int> (3);
// auto b = vp.give<double>(5, 7);
// auto c = vp.give<std::tuple<short, short, short> >(2);
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
// auto v = vp.give< vec<vec<int>> > (3); // A vector of vector of vectors.
// for (int i = 0; i < 3; ++i)
// {
//   vp.give<vec<int>>(i * 2 + 1).swap(v[i]);
//   for (int j = 0, jend = v[i].size(); j < jend; ++j)
//     vp.give<int>(j + 5).swap(v[i][j]);
// }




}