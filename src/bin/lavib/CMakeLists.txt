# ================================================================
# bfsca - lavib CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(BIN_NAME lavib.x)

set(SRC_LIST 
  lavib.f  rclosure.f
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
target_link_libraries(${BIN_NAME} nr)
target_link_libraries(${BIN_NAME} quadpack)

# ================================================================
# Installation
install(TARGETS ${BIN_NAME} DESTINATION ${INSTA_BIN})

# ================================================================
