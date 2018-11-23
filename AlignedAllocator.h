#pragma once

#include <mm_malloc.h>
#include <vector>

////////////////////////////////////////////////////////////////////
// A lattice of something, but assume the something is SIMDized.
////////////////////////////////////////////////////////////////////
#define GRID_ALLOC_ALIGN (2*1024*1024)

template<typename _Tp>
class alignedAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef alignedAllocator<_Tp1> other; };
  alignedAllocator() throw() { }
  alignedAllocator(const alignedAllocator&) throw() { }
  template<typename _Tp1> alignedAllocator(const alignedAllocator<_Tp1>&) throw() { }
  ~alignedAllocator() throw() { }
  pointer       address(reference __x)       const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
    size_type bytes = __n*sizeof(_Tp);
    pointer ptr = (_Tp *) _mm_malloc(bytes,GRID_ALLOC_ALIGN);
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) { 
    _mm_free((void *)__p); 
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return false; }
template<class T> using Vector     = std::vector<T,alignedAllocator<T> >;           
    
