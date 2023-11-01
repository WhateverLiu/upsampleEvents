namespace Charlie {


// =============================================================================
// Swap two std::vectors of different types. After swapping, both vectors'
//   sizes will become 0. For each of them, the capacity pointer could be shifted
//   to ensure that the capacity is a multiple of sizeof(type). 
//   No allocation will occur. The vectors maintain their original types.
// The function is thread-safe as it has no data member.
// =============================================================================
struct DarkSwapVec
{
  DarkSwapVec()
  {
    static_assert( sizeof(std::size_t) == sizeof(void*) );
    
    
    // Test if std::vector is implemented in the way that can be exploited.
    constexpr const unsigned fullSize = 31;
    constexpr const unsigned subSize  = 13;
    
    
    typedef std::tuple<char, char, char> tupe; // Arbitrarily selected type.
    std::vector<tupe> a(fullSize); // Just for creating a vector of an arbitrary type.
    a.resize(subSize); // Downsize the vector.
    
    
    // Vector header must contain exactly 3 pointers, e.g. p1, p2, p3.
    static_assert( sizeof(a) == sizeof(std::size_t) * 3 );
    auto ptr = (std::size_t*)(&a); // Read but won't write.
    bool sizePointerCool     = ptr[0] + sizeof(tupe) * subSize  == ptr[1];
    bool capacityPointerCool = ptr[0] + sizeof(tupe) * fullSize == ptr[2];
    
    
    if (!sizePointerCool or !capacityPointerCool)
      throw std::runtime_error("Layout of std::vector header is unsupported.");
  } 
  
  
  template<typename S, typename T>
  void operator()(std::vector<S> &x, std::vector<T> &y)
  { 
    // Enforce destructor on elements if type is nontrivial.
    //   Will not trigger reallocation.
    x.clear();
    y.clear();
    if constexpr ( std::is_same<S,T>::value ) 
    {
      x.swap(y);
    }
    else
    {
      static_assert( sizeof(std::size_t) == sizeof(void*) );
      static_assert( sizeof(x) == sizeof(y) );
      static_assert( sizeof(x) == sizeof(std::size_t) * 3 );
      std::size_t xt[3], yt[3];
      mmcp(xt, &x, sizeof(x));
      mmcp(yt, &y, sizeof(y));
      xt[2] -= (xt[2] - xt[0]) % sizeof(T);
      yt[2] -= (yt[2] - yt[0]) % sizeof(S);
      mmcp(&x, yt, sizeof(x));
      mmcp(&y, xt, sizeof(y));
      // Do not overthink. Optmization over whether S or T is byte is not worth it.  
    }
  }
};



}
