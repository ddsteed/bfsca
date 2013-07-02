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
