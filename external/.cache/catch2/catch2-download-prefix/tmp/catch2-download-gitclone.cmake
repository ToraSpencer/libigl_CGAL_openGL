# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitclone-lastrun.txt" AND EXISTS "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitinfo.txt" AND
  "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitclone-lastrun.txt" IS_NEWER_THAN "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/catch2"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/catch2'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "D:/Git/cmd/git.exe"
            clone --no-checkout --config "advice.detachedHead=false" --config "advice.detachedHead=false" -c http.sslVerify=false "https://github.com/catchorg/Catch2.git" "catch2"
    WORKING_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/catchorg/Catch2.git'")
endif()

execute_process(
  COMMAND "D:/Git/cmd/git.exe"
          checkout "v2.11.0" --
  WORKING_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/catch2"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'v2.11.0'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "D:/Git/cmd/git.exe" -c;http.sslVerify=false
            submodule update --recursive --init 
    WORKING_DIRECTORY "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/catch2"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/catch2'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitinfo.txt" "G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'G:/gitRepositories/libigl_CGAL_openGL/external/.cache/catch2/catch2-download-prefix/src/catch2-download-stamp/catch2-download-gitclone-lastrun.txt'")
endif()
