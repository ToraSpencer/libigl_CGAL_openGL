cmake_minimum_required(VERSION 3.12)
 

# 0. blah blah blah
project(myTest VERSION 1.0)                 # 设定项目名称是helloWorld，版本是1.0
set(CMAKE_CXX_STANDARD 11)                          # 指定项目的C++标准：
set(CMAKE_CXX_STANDARD_REQUIRED True)


message(STATUS "CMAKE_SYSTEM ==  ${CMAKE_SYSTEM}")
message(STATUS "CMAKE_SYSTEM_NAME ==  ${CMAKE_SYSTEM_NAME}")


# 1. 导入源代码文件：

#       aux_source_directory()——查找某路径下所有源文件：
aux_source_directory(${PROJECT_SOURCE_DIR} SRC_FILES)    # aux_source_directory(路径，存储源文件列表的变量)



# 2. 手动链接静态库：

#   2.1 包含静态库头文件路径
include_directories(${PROJECT_SOURCE_DIR}/../include)
include_directories(${PROJECT_SOURCE_DIR}/../external)
include_directories(${PROJECT_SOURCE_DIR}/../external/eigen)
include_directories(${PROJECT_SOURCE_DIR}/../external/glad/include)
include_directories(${PROJECT_SOURCE_DIR}/../external/glfw/include)


#   2.2 link_directories()——包含需要链接的库文件路径
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")

    # 设定只生成Release版本：
    set(CMAKE_CONFIGURATION_TYPES "Release")        # 配置版本；默认情形下为Debug;Release;MinSizeRel;RelWithDebInfo
    message(STATUS "CMAKE_CONFIGURATION_TYPES == ${CMAKE_CONFIGURATION_TYPES}")

    link_directories(${PROJECT_SOURCE_DIR}/lib/WINDOWS/${CMAKE_CONFIGURATION_TYPES})

ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Apple")

ENDIF (CMAKE_SYSTEM_NAME MATCHES "Windows")
 

# 3. 预处理器

#   add_definitions(-DXXX)——预处理器中添加宏XXX
add_definitions(-D_CRT_SECURE_NO_WARNINGS)         


# 4. 设定输出内容：

#   4.1 设定输出的可执行文件名称为"myTest"; 对应的源代码列表是SRC_FILES中存储的值；
add_executable(myTest ${SRC_FILES})

#   4.2 target_link_libraries()——用在add_executable()之后
target_link_libraries(myTest  igl)
target_link_libraries(myTest  glad)
target_link_libraries(myTest  glfw3)
target_link_libraries(myTest  igl_opengl)
target_link_libraries(myTest  igl_opengl_glfw)