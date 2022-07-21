# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 3.1)

project(gmp-download NONE)

include(ExternalProject)
ExternalProject_Add(gmp-download
                    TLS_VERIFY;OFF;GIT_CONFIG;advice.detachedHead=false;URL;https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/GMP/5.0.1/gmp-all-CGAL-3.9.zip;URL_MD5;508c1292319c832609329116a8234c9f
                    SOURCE_DIR          "G:/gitRepositories/libigl_CGAL_openGL/cmake/../external/gmp"
                    BINARY_DIR          "G:/gitRepositories/libigl_CGAL_openGL/build/gmp-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
