#ifndef IGL_DECIMATE_TRIVIAL_CALLBACKS_H
#define IGL_DECIMATE_TRIVIAL_CALLBACKS_H
#include "igl_inline.h"
#include "decimate_callback_types.h"


namespace igl
{

  /*
      Function to build trivial preand post collapse actions.
  
       Outputs:
         always_try       function that always returns true (always attempt the next edge collapse)
         never_care      fuction that is always a no-op (never have a post collapse response)
    */
  IGL_INLINE void decimate_trivial_callbacks(
    decimate_pre_collapse_callback  & always_try,
    decimate_post_collapse_callback & never_care);
};


#ifndef IGL_STATIC_LIBRARY
#  include "decimate_trivial_callbacks.cpp"
#endif

#endif 

