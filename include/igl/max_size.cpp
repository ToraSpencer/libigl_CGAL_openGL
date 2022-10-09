// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "max_size.h"


template <typename T>
IGL_INLINE int igl::max_size(const std::vector<T> & vers)
{
  int max_size = -1;
  for(
    typename std::vector<T>::const_iterator iter = vers.begin();
    iter != vers.end(); 
    iter++)
  {
    int size = (int)iter->size();
    max_size = (max_size > size ? max_size : size);
  }
  return max_size;
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template int igl::max_size<std::vector<long, std::allocator<long> > >(std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > > const&);
// generated by autoexplicit.sh
template int igl::max_size<std::vector<unsigned int, std::allocator<unsigned int> > >(std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > > const&);
// generated by autoexplicit.sh
template int igl::max_size<std::vector<int, std::allocator<int> > >(std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&);
// generated by autoexplicit.sh
template int igl::max_size<std::vector<double, std::allocator<double> > >(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&);
template int igl::max_size<std::vector<bool, std::allocator<bool> > >(std::vector<std::vector<bool, std::allocator<bool> >, std::allocator<std::vector<bool, std::allocator<bool> > > > const&);
template int igl::max_size<std::vector<float, std::allocator<float> > >(std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > const&);
#ifdef WIN32
template int igl::max_size<class std::vector<__int64,class std::allocator<__int64> > >(class std::vector<class std::vector<__int64,class std::allocator<__int64> >,class std::allocator<class std::vector<__int64,class std::allocator<__int64> > > > const &);
#endif
#endif
