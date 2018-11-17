cmake_minimum_required(VERSION 3.1)



################################################################################

### Configuration

set(LIBHEDRA_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(LIBHEDRA_SOURCE_DIR "${LIBHEDRA_ROOT}/include")
set(LIBHEDRA_EXTERNAL "${LIBHEDRA_ROOT}/external")

set(LIBIGL_ROOT "${LIBHEDRA_EXTERNAL}/libigl")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

include_directories(${LIBHEDRA_SOURCE_DIR})
