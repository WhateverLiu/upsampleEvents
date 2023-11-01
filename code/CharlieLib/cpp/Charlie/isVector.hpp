#pragma once


namespace Charlie {

// =============================================================================
// Compile time computation to decide if type is a vector of any type.
// =============================================================================
// template<typename T>
// struct isVector: public std::false_type {};
// 
// template<typename T, typename A>
// struct isVector<std::vector<T, A> >: public std::true_type {};
// // Usage: if constexpr (isVector<DataType>{}) do something.


template<typename T>
struct isVector { constexpr bool operator()() { return false; } };


template<typename T, typename A> // Partial template specialization.
struct isVector<std::vector<T, A>> 
{
  constexpr bool operator()() { return true; };
};
// =============================================================================



}

