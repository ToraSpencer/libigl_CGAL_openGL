# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 3.1)

project(test_data-download NONE)

include(ExternalProject)
ExternalProject_Add(test_data-download
                    TLS_VERIFY;OFF;GIT_CONFIG;advice.detachedHead=false;GIT_REPOSITORY;https://github.com/libigl/libigl-tests-data;GIT_TAG;19cedf96d70702d8b3a83eb27934780c542356fe
                    SOURCE_DIR          "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/../tests/data"
                    BINARY_DIR          "G:/gitRepositories/libigl_CGAL_openGL/build/test_data-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
