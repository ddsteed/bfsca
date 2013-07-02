# ================================================================
# bfsca - src/lib/libqt CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(LIB_NAME qt)

set(LIBMINTS_SRC
  clebg.f  ylm.f  pml.f
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
