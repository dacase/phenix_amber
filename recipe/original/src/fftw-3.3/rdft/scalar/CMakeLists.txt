set(RDFT_SCALAR_SOURCES hfb.c r2c.c r2r.c hc2c.c)

include_directories(.)

add_library(rdft_scalar OBJECT ${RDFT_SCALAR_SOURCES})
make_pic_if_needed(rdft_scalar)

add_subdirectory(r2cb)
add_subdirectory(r2cf)
add_subdirectory(r2r)
