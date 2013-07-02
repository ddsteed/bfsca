# ================================================================
# bfsca - sumxsc CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(BIN_NAME sumxsc.x)

set(SRC_LIST 
  sumxsc.f
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
