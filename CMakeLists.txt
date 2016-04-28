cmake_minimum_required(VERSION 2.8)
project(convex_decomposer)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set (HACD_PATH "" CACHE PATH "Path to HACD")
set (HACD_PATH "" CACHE PATH "Path to HACD")

find_package(Bullet REQUIRED)

if(BULLET_FOUND)

find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)

include_directories(include ${Eigen_INCLUDE_DIRS} ${BULLET_INCLUDE_DIR} ${HACD_PATH})

add_executable(convex_decomposer
                  src/convex_decomposer.cpp
                  src/cd_wavefront.cpp
                  ${HACD_PATH}/hacdGraph.cpp
                  ${HACD_PATH}/hacdHACD.cpp
                  ${HACD_PATH}/hacdICHull.cpp
                  ${HACD_PATH}/hacdManifoldMesh.cpp)

target_link_libraries(convex_decomposer ${BULLET_LIBRARIES})

else(BULLET_FOUND)

message(STATUS "Bullet has not been found! Skipping the dynamics simulator...")

endif()
