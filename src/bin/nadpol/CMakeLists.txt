# ================================================================
# bfsca - nadpol CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(BIN_NAME nadpol.x)

set(SRC_LIST 
  nadpol.f
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
target_link_libraries(${BIN_NAME} qt)

# ================================================================
# Installation
install(TARGETS ${BIN_NAME} DESTINATION ${INSTA_BIN})

# ================================================================
