#pragma once


namespace Charlie
{
template <int size> // Number of items to print, e.g. 100.
struct ProgressBar // Suitable in multithreading environment.
{
  int64_t which;
  std::size_t gap;
  bool p[size];
  ProgressBar(std::size_t imax)
  { 
    std::fill(p, p + size, false); 
    gap = (imax + (size - 1)) / size;
  }
  int64_t operator()(std::size_t i)
  {
    which = i / gap;
    if (!p[which]) 
    {
      p[which] = true;
      return which;
    }
    return -1;
  }
};
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

