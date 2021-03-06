# Install script for directory: libigl_CGAL_openGL/tutorial

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/libigl_tutorials")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("./glad/cmake_install.cmake")
  include("./glfw/cmake_install.cmake")
  include("./101_FileIO/cmake_install.cmake")
  include("./102_DrawMesh/cmake_install.cmake")
  include("./103_Events/cmake_install.cmake")
  include("./104_Colors/cmake_install.cmake")
  include("./107_MultipleMeshes/cmake_install.cmake")
  include("./108_MultipleViews/cmake_install.cmake")
  include("./110_MshView/cmake_install.cmake")
  include("./201_Normals/cmake_install.cmake")
  include("./202_GaussianCurvature/cmake_install.cmake")
  include("./203_CurvatureDirections/cmake_install.cmake")
  include("./204_Gradient/cmake_install.cmake")
  include("./205_Laplacian/cmake_install.cmake")
  include("./206_GeodesicDistance/cmake_install.cmake")
  include("./207_PolygonLaplacian/cmake_install.cmake")
  include("./301_Slice/cmake_install.cmake")
  include("./302_Sort/cmake_install.cmake")
  include("./303_LaplaceEquation/cmake_install.cmake")
  include("./304_LinearEqualityConstraints/cmake_install.cmake")
  include("./305_QuadraticProgramming/cmake_install.cmake")
  include("./306_EigenDecomposition/cmake_install.cmake")
  include("./401_BiharmonicDeformation/cmake_install.cmake")
  include("./402_PolyharmonicDeformation/cmake_install.cmake")
  include("./403_BoundedBiharmonicWeights/cmake_install.cmake")
  include("./404_DualQuaternionSkinning/cmake_install.cmake")
  include("./405_AsRigidAsPossible/cmake_install.cmake")
  include("./406_FastAutomaticSkinningTransformations/cmake_install.cmake")
  include("./407_BiharmonicCoordinates/cmake_install.cmake")
  include("./408_DirectDeltaMush/cmake_install.cmake")
  include("./501_HarmonicParam/cmake_install.cmake")
  include("./502_LSCMParam/cmake_install.cmake")
  include("./503_ARAPParam/cmake_install.cmake")
  include("./507_Planarization/cmake_install.cmake")
  include("./701_Statistics/cmake_install.cmake")
  include("./702_WindingNumber/cmake_install.cmake")
  include("./703_Decimation/cmake_install.cmake")
  include("./704_SignedDistance/cmake_install.cmake")
  include("./705_MarchingCubes/cmake_install.cmake")
  include("./707_SweptVolume/cmake_install.cmake")
  include("./708_Picking/cmake_install.cmake")
  include("./709_SLIM/cmake_install.cmake")
  include("./711_Subdivision/cmake_install.cmake")
  include("./712_DataSmoothing/cmake_install.cmake")
  include("./713_ShapeUp/cmake_install.cmake")
  include("./715_MeshImplicitFunction/cmake_install.cmake")
  include("./716_HeatGeodesics/cmake_install.cmake")
  include("./718_IterativeClosestPoint/cmake_install.cmake")
  include("./719_ExplodedView/cmake_install.cmake")
  include("./720_BlueNoise/cmake_install.cmake")
  include("./721_VectorFieldSmoothing/cmake_install.cmake")
  include("./722_VectorParallelTransport/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "./${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
