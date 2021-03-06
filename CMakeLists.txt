cmake_minimum_required (VERSION 3.5)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(Embedding)

set(CMAKE_CXX_STANDARD 14)

#Output directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

#grab important packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


#Gauss
set(LIBIGL_INCLUDE_PATH "/home/vismay/recode/libigl/include") # Weird that I have to do this for GAUSS seperately
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fPIC)
set(Gauss_ROOT_DIR "/home/vismay/recode/GAUSS/")
find_package(GAUSS REQUIRED)

#Libigl
set(LIBIGL_INCLUDE_DIR "/home/vismay/recode/libigl/include")
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           ON)
option(LIBIGL_WITH_PNG              "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML              "Use XML"            OFF)

find_package(LIBIGL REQUIRED QUIET)



#get the source code files
file(GLOB SRC_CODE ${Embedding_SOURCE_DIR}/src/*.cpp)

include_directories(${Gauss_INCLUDE_DIRS})

#this call is necessary if any targets are GUI applications
#initUI()  #Is this important????

#Example codes
link_directories(${Gauss_LIB_DIR_RELEASE})
add_executable(Embedding ${Embedding_SOURCE_DIR}/src/exampleK.cpp)
target_link_libraries(Embedding ${Gauss_LIBS} ${UI_LIBS} -lpthread igl::core igl::opengl igl::opengl_glfw)
