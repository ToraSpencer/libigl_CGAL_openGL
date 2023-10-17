# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/triangle"
  "G:/gitRepositories/libigl_CGAL_openGL/build32/triangle-build"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/tmp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/src/triangle-download-stamp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/src"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/src/triangle-download-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/src/triangle-download-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/triangle/triangle-download-prefix/src/triangle-download-stamp${cfgdir}") # cfgdir has leading slash
endif()