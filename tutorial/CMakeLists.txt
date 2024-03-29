cmake_minimum_required(VERSION 3.16)
project(libhedra_tutorials)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### Choose which chapters to compile
option(TUTORIALS_CHAPTER1 "Compile chapter 1" ON)
option(TUTORIALS_CHAPTER2 "Compile chapter 2" ON)
option(TUTORIALS_CHAPTER3 "Compile chapter 3" ON)
option(TUTORIALS_CHAPTER4 "Compile chapter 4" ON)

### libIGL options:
option(LIBIGL_EMBREE           "Build target igl::embree"           ON)
option(LIBIGL_GLFW             "Build target igl::glfw"             ON)
option(LIBIGL_IMGUI            "Build target igl::imgui"            ON)
option(LIBIGL_OPENGL           "Build target igl::opengl"           ON)
option(LIBIGL_PNG              "Build target igl::png"              ON)
### set(LIBIGL_WITH_PNG CACHE BOOL ON)


### Adding libIGL and libhedra: choose the path to your local copy
include(libigl)
include(libhedra)

if (LIBHEDRA_WITH_CERES)
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
endif()

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()




# Store location of tutorial/shared directory
set(TUTORIAL_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shared CACHE PATH "location of shared tutorial resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})



# Chapter 1
if(TUTORIALS_CHAPTER1)
  add_subdirectory("101_Visualization")
  add_subdirectory("102_Evaluation")
endif()

# Chapter 2
if(TUTORIALS_CHAPTER2 AND LIBHEDRA_WITH_CERES)
  add_subdirectory("201_ModelingAffine")
  add_subdirectory("202_PolyhedralPatterns")
endif()

# Chapter 3
if(TUTORIALS_CHAPTER3 AND LIBHEDRA_WITH_CERES)
  add_subdirectory("301_ComplexMoebiusDeformation")
  add_subdirectory("302_ComplexMoebiusInterpolation")
  add_subdirectory("303_QuatMoebiusDeformation")
  add_subdirectory("304_RegularMeshes")
endif()

# Chapter 4
if(TUTORIALS_CHAPTER4)
  add_subdirectory("401_PolygonalSubdivision")
  add_subdirectory("402_CanonicalMoebiusSubdivision")
endif()


