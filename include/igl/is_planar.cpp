// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "is_planar.h"
IGL_INLINE bool igl::is_planar(const Eigen::MatrixXd & vers)
{
  if(vers.size() == 0) return false;
  if(vers.cols() == 2) return true;
  for(int i = 0;i<vers.rows();i++)
  {
    if(vers(i,2) != 0) return false;
  }
  return true;
}
