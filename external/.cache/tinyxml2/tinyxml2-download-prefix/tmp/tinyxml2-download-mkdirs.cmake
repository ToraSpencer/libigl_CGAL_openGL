# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/tinyxml2"
  "G:/gitRepositories/libigl_CGAL_openGL/build32/tinyxml2-build"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/tmp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/src/tinyxml2-download-stamp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/src"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/src/tinyxml2-download-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/src/tinyxml2-download-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/tinyxml2/tinyxml2-download-prefix/src/tinyxml2-download-stamp${cfgdir}") # cfgdir has leading slash
endif()
