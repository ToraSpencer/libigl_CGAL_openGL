## Copyright 2009-2020 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

# use default install config
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/embree-config-install.cmake)

# and override path variables to match for build directory
SET(EMBREE_INCLUDE_DIRS G:/gitRepositories/libigl_CGAL_openGL/external/embree/include)
SET(EMBREE_LIBRARY G:/gitRepositories/libigl_CGAL_openGL/build32/embree/embree3.lib)
SET(EMBREE_LIBRARIES ${EMBREE_LIBRARY})
