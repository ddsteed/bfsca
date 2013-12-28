# ================================================================
# bfsca - src/lib/libquadpack CMakeLists.txt
#
# Hao Feng (RDS), XHU
#
# ================================================================
set(LIB_NAME quadpack)

set(LIBMINTS_SRC
  d1mach.f	dqage.f		dqagse.f	dqk15.f		
  dqk31.f	dqk51.f		dqpsrt.f    dqag.f	        
  dqags.f	dqelg.f		dqk21.f		dqk41.f		
  dqk61.f
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
