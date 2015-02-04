# ============================================================
# BFSCA CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ============================================================
# Primary Make
cmake_minimum_required(VERSION 2.8)

# Project name and support languages
project(BFSCA
        CXX Fortran C)

# handle mem > 2 GB by using Intel compiler
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS "-shared-intel -mcmodel medium")
   endif()
endif()

# Source code subdirectory
add_subdirectory(src)

# ============================================================
# Package
set(CPACK_GENERATOR "TGZ") 

# Setup package version
set(CPACK_PACKAGE_VERSION_MAJOR 0) 
set(CPACK_PACKAGE_VERSION_MINOR 1) 
set(CPACK_PACKAGE_VERSION_PATCH 0)

# ’call’ CPack
include(CPack)

# ============================================================
