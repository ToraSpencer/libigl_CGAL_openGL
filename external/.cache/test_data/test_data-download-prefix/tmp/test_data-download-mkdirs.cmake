# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/../tests/data"
  "G:/gitRepositories/libigl_CGAL_openGL/build32/test_data-build"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/tmp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/src/test_data-download-stamp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/src"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/src/test_data-download-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/src/test_data-download-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/test_data/test_data-download-prefix/src/test_data-download-stamp${cfgdir}") # cfgdir has leading slash
endif()
