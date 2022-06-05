
if(NOT "E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitinfo.txt" IS_NEWER_THAN "E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "E:/gitRepositories/libigl_CGAL_openGL/cmake/../external/cgal"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'E:/gitRepositories/libigl_CGAL_openGL/cmake/../external/cgal'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "D:/Git/cmd/git.exe" -c http.sslVerify=false clone --no-checkout --config "advice.detachedHead=false" --config "advice.detachedHead=false" "https://github.com/CGAL/cgal.git" "cgal"
    WORKING_DIRECTORY "E:/gitRepositories/libigl_CGAL_openGL/cmake/../external"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/CGAL/cgal.git'")
endif()

execute_process(
  COMMAND "D:/Git/cmd/git.exe" -c http.sslVerify=false checkout f7c3c8212b56c0d6dae63787efc99093f4383415 --
  WORKING_DIRECTORY "E:/gitRepositories/libigl_CGAL_openGL/cmake/../external/cgal"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'f7c3c8212b56c0d6dae63787efc99093f4383415'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "D:/Git/cmd/git.exe" -c http.sslVerify=false submodule update --recursive --init 
    WORKING_DIRECTORY "E:/gitRepositories/libigl_CGAL_openGL/cmake/../external/cgal"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'E:/gitRepositories/libigl_CGAL_openGL/cmake/../external/cgal'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitinfo.txt"
    "E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'E:/gitRepositories/libigl_CGAL_openGL/external/.cache/cgal/cgal-download-prefix/src/cgal-download-stamp/cgal-download-gitclone-lastrun.txt'")
endif()

