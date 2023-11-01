

namespace Charlie {

#define XXH_STATIC_LINKING_ONLY   /* access advanced declarations */
#define XXH_IMPLEMENTATION   /* access definitions */
#include "xxhash.h"


// XXH64(&v[0], sizeof(int) * v.size(), 42);
inline std::size_t hash(const void *x, std::size_t Nbyte, std::size_t seed)
{
  return XXH64(x, Nbyte, seed);
}


#undef XXH_IMPLEMENTATION
#undef XXH_STATIC_LINKING_ONLY




}