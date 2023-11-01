namespace Charlie {
// =============================================================================
// Best practice for using temporary vector allocator:
//
//   NEVER DO ANY OPERATION TAHT WILL EXPAND THE CAPACITY OF A VECTOR THAT
//   ALREADY HAS NONZERO CAPACITY. HERE CAPACITY = std::vector::capacity().
//
// Violation to the above recommendation will not lead to errors, but could
//   substantially fragment the memory pool.
// 
// For example, for finding the positive elements in vector `x`, one should much
//   prefer:
// 
//   std::vector<T, A> xPositive(x.size()); 
//   xPositive.resize(0);
//   for (int i = 0, iend = x.size(); i < iend; ++i)
//     if (x[i] > 0) xPositive.emplace_back(x[i]);
//  
// than naively:
// 
//   std::vector<T, A> xPositive;
//   for (int i = 0, iend = x.size(); i < iend; ++i)
//     if (x[i] > 0) xPositive.emplace_back(x[i]);
// =============================================================================
#define vec std::vector


/*
template <typename T>
struct RawVec
{
  T *begin, *end, *cap;
  RawVec(){ begin = end = cap = nullptr; }
  void reset(std::size_t size, std::size_t alignment) // alignment MUST be power of 2.
  {
    begin = (T*)std::malloc(size * sizeof(T));
    cap = end = begin + size; 
  }
  RawVec(std::size_t size) { reset(size); }
  ~RawVec() 
  { 
    for (end -= 1; end >= begin; --end) end->~T();
    free(begin); 
    begin = end = cap = nullptr; 
  }
  void resize(std::size_t size, T &&t = T())
  {
    std::size_t oldSize = end - begin;
    for (auto it = begin + size; it < end; ++it) it->~T();
    
    if (size == oldSize) return;
    if (size < oldSize)
    {
      for (auto it = begin + size; it != end; ++it)
        
    }
    if (size <= cap - begin) end = begin + size;
    else { free(begin); reset(size); }
  }
  void emplace_back(T &&x)
  {
    
  }
};
*/


template <std::size_t alignment = 16>
struct SAAmemPool // Single threaded
{
  std::size_t initialMemBlockByte;
  vec<char> *Xback;
  vec<vec<char>> X;
  
  
  std::size_t byteAllocated()
  {
    std::size_t S = 0;
    for (auto &x: X) S += x.capacity();
    return S;
  }
  
  
  // Find the next char address that is 16-byte aligned.
  char *nextStart (char *x)
  {
    constexpr const std::size_t m = alignment - 1;
    return (char*)((std::size_t(x) + 8 + m) & m);
  }
  
  
  bool isAligned(char *x) { return (std::size_t(x) & (alignment - 1)) == 0; }
  
  
  void reset(std::size_t&& initialMemBlockByte)
  {
    initialMemBlockByte = std::max(initialMemBlockByte, sizeof(std::size_t) * 8);
    this->initialMemBlockByte = initialMemBlockByte;
    X.resize(0);
    X.emplace_back(vec<char> (initialMemBlockByte));
    Xback = &X.back();
    if (isAligned(Xback->data())) std::runtime_error("Storage start is not " + 
        std::to_string(alignment) + "-byte aligned.\n");
  }
  
  
  SAAmemPool() { reset(32 * 1024 * 1024); }
  SAAmemPool(std::size_t&& initialMemBlockByte) { reset(initialMemBlockByte); }
  
  
  // &*x.end() is always a multiple of 16 bytes.
  char *allocate(std::size_t size) // Always return 16-byte aligned address.
  {
    if (size == 0) return nullptr;
    auto &x = *Xback;
    auto nst = &*x.end(); // New begin.
    auto nend = nextStart(nst + size); // New end, also 16-byte aligned.
    if (nend < x.data() + x.capacity())
    {
      x.resize(nend - x.data());
      std::memset(nend - sizeof(std::size_t), 0, sizeof(std::size_t)); // This block is not registered for deletion.
      return nst;
    }
    X.emplace_back(vec<char> (
        std::max(X.back().capacity() * 2, size + initialMemBlockByte)
    ));
    Xback = &X.back();
    if (isAligned(Xback->data())) std::runtime_error("Storage start is not " + 
        std::to_string(alignment) + "-byte aligned.\n");
    return allocate(size);
  }
  
  
  std::size_t regis;
  void autoDeallocate()
  {
    auto &x = *Xback;
    if ( x.size() )
    {
      std::memcpy(&regis, &*x.end() - sizeof(regis), sizeof(regis));
      if (regis == 0) return;
      x.resize(x.size() - regis);
      autoDeallocate();
    }
    else
    {
      if (Xback <= X.data()) return;
      Xback -= 1;
      autoDeallocate();
    }
  }
  
  
  // ptr is guaranteed to be 16-byte aligned.
  void deallocate(char *ptr, std::size_t size) noexcept
  {
    char *p = nextStart(ptr + size); // End of the block.
    if (p == &*Xback->end()) // It IS the last block to be released.
    {
      Xback->resize(ptr - Xback->data());
      autoDeallocate();
    }
    else
    {
      // Write the block's byte size to the block's last 8 bytes.
      regis = p - ptr;
      std::memcpy(p - sizeof(regis), &regis, sizeof(regis));
    }
  }
  
  
};


/*
template <std::size_t alignment = 16>
struct SAAmemPool
{
  vec<SAAmemPoolST> P;
  SAAmemPool(){}
  SAAmemPool(std::size_t&& initialMemBlockByte)
  {
  }
};
*/




// Stack alike allocator
template <typename T>
struct SAA
{
  SAAmemPool<> *pool;
  SAA() noexcept { pool = nullptr; }
  SAA(SAAmemPool<> *pool): pool(pool) {}
  typedef T value_type;
  template<class U> SAA(const SAA<U>&) noexcept {}
  template<class U> bool operator==(const SAA<U>&) const noexcept { return true; }
  template<class U> bool operator!=(const SAA<U>&) const noexcept { return false; }
  T *allocate(std::size_t n) { return reinterpret_cast<T*> (pool->allocate(n)); }
  void deallocate(T *ptr, std::size_t n) noexcept
  {
    pool->deallocate(reinterpret_cast<char*>(ptr), n * sizeof(T));
  }
};




// =============================================================================
// Temporary vector generation.
// =============================================================================
template<typename T> using mvec = vec<T, SAA<T> >;


template <typename T>
mvec<T> mvecMakeCore(std::size_t size, SAAmemPool<> &pool)
{
  SAA<T> A(&pool);
  return mvec<T>(size, A);
}


// Undefined behavior if pool is nullptr.
template <typename T, typename... Args>
auto mvecMake(SAAmemPool<> &pool, std::size_t size, Args... restSizes)
{
  if constexpr (sizeof...(restSizes) == 0)
    return mvecMakeCore<T>(size, pool);
  else
  {
    auto rst = mvecMake<decltype( mvecMake<T>(pool, restSizes...) )>(pool, size);
    for (auto &x: rst) mvecMake<T>(pool, restSizes...).swap(x);
    return rst;
  }
}






#undef vec























}