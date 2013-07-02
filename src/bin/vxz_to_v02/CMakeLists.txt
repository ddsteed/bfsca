# ================================================================
# bfsca - vxz_to_v02 CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(BIN_NAME vxz_to_v02.x)

set(SRC_LIST 
  vxz_to_v02.f
  )

# ================================================================
# executable files
add_executable(${BIN_NAME}  ${SRC_LIST})

# ================================================================
# include directory
#include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/src/lib)

# ================================================================
# link library
#target_link_libraries(${BIN_NAME} BFSCA_util)

# ================================================================
# Installation
install(TARGETS ${BIN_NAME} DESTINATION ${INSTA_BIN})

# ================================================================
