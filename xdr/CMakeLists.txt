# Policy
cmake_minimum_required( VERSION 3.0) 



if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()


include(ConfigureChecks.cmake)
add_definitions(-DHAVE_CONFIG_H=1)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

include_directories("${PROJECT_BINARY_DIR}")
add_subdirectory(src)
