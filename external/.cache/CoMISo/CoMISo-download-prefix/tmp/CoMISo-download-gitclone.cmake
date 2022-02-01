
if(NOT "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitinfo.txt" IS_NEWER_THAN "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external/CoMISo"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external/CoMISo'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/local/bin/git" -c http.sslVerify=false clone --no-checkout --config "advice.detachedHead=false" --config "advice.detachedHead=false" "https://github.com/libigl/CoMISo.git" "CoMISo"
    WORKING_DIRECTORY "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/libigl/CoMISo.git'")
endif()

execute_process(
  COMMAND "/usr/local/bin/git" -c http.sslVerify=false checkout d60aa4759fba76b0b793b1efb090b7a771dd7c56 --
  WORKING_DIRECTORY "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external/CoMISo"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'd60aa4759fba76b0b793b1efb090b7a771dd7c56'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/local/bin/git" -c http.sslVerify=false submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external/CoMISo"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/wuhan/gitRepositories/libigl_CGAL_openGL/cmake/../external/CoMISo'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitinfo.txt"
    "/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/wuhan/gitRepositories/libigl_CGAL_openGL/external/.cache/CoMISo/CoMISo-download-prefix/src/CoMISo-download-stamp/CoMISo-download-gitclone-lastrun.txt'")
endif()

