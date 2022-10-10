#ifndef IGL_MAX_FACES_STOPPING_CONDITION_H
#define IGL_MAX_FACES_STOPPING_CONDITION_H
#include "igl_inline.h"
#include "decimate_callback_types.h"
#include <Eigen/Core>
#include <vector>
#include <set>
#include <functional>

namespace igl
{
    // 生成decimate()接口中迭代终止函数子stopping_condition();
    /*
       Stopping condition function compatible with igl::decimate. The outpute
       function handle will return true if number of faces is less than max_m
  
       Inputs:
         m                   reference to working variable initially should be set to current number of faces.
         orig_m            number (size) of original face list _**not**_ including any
                                       faces added to handle phony boundary faces connecting to point at
                                       infinity. For closed meshes it's safe to set this to tris.rows()
         max_m           maximum number of faces

       Outputs:
         stopping_condition
  */
  IGL_INLINE void max_faces_stopping_condition(
    int & m,
    const int orig_m,
    const int max_m,
    decimate_stopping_condition_callback & stopping_condition);
  IGL_INLINE decimate_stopping_condition_callback
    max_faces_stopping_condition(
      int & m,
      const int orign_m,
      const int max_m);
}

#ifndef IGL_STATIC_LIBRARY
#  include "max_faces_stopping_condition.cpp"
#endif
#endif

