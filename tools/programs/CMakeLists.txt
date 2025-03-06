#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.15)



PROJECT(dBased)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON) 

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(dBased_CMAKE_MODULES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
  add_definitions(/MP)
endif()

find_package(OpenMP)

set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_MINSIZEREL_POSTFIX "s")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")

include_directories(${DISCAMB_INCLUDE_PATH})
include_directories(${NLOHMANN_JSON_PATH})
include_directories(${EIGEN_PATH})
link_directories(${DISCAMB_LIB_PATH})

if(WIN32)
    SET(CMAKE_DETECTED_WIN32 TRUE)
endif(WIN32)

ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/src")

set (DISCAMBBASED_BUILD_DIR "${PROJECT_BINARY_DIR}/build")


