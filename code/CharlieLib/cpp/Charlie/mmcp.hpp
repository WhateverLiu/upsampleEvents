#pragma once


namespace Charlie {

template<typename DestPtr, typename SourPtr>
inline void mmcp(DestPtr dest, SourPtr src, std::size_t Nbyte)
{ 
  std::memcpy((void*)dest, (void*)src, Nbyte);
}


}