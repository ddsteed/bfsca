# ================================================================
# bfsca - src/lib/libnr CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(LIB_NAME nr)

set(LIBMINTS_SRC
  spline.f   splint.f  gaussian.f  gauleg.f
)

# ================================================================
# Library files
add_library(${LIB_NAME} STATIC ${LIBMINTS_SRC})

# ================================================================
# include directory
#include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/src/lib)

# ================================================================
# Installation
install(TARGETS ${LIB_NAME} ${LIB_NAME}
        LIBRARY DESTINATION ${INSTA_LIB}
        ARCHIVE DESTINATION ${INSTA_LIB})
# ================================================================
