# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/libigl-imgui"
  "G:/gitRepositories/libigl_CGAL_openGL/build32/libigl-imgui-build"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/tmp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/src/libigl-imgui-download-stamp"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/src"
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/src/libigl-imgui-download-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/src/libigl-imgui-download-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/libigl-imgui/libigl-imgui-download-prefix/src/libigl-imgui-download-stamp${cfgdir}") # cfgdir has leading slash
endif()
