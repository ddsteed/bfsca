# ================================================================
# bfsca src CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
add_subdirectory(lib)
add_subdirectory(bin)

# ================================================================
# add openmp support
# FIND_PACKAGE(OpenMP)
# IF(OPENMP_FOUND)
#   SET(CMAKE_C_FLAGS           "${CMAKE_C_FLAGS}          ${OpenMP_C_FLAGS}")
#   SET(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
#   SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# ENDIF()

# ================================================================
