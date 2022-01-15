include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(dft_scalar OBJECT n.c t.c)
make_pic_if_needed(dft_scalar)

add_subdirectory(codelets)