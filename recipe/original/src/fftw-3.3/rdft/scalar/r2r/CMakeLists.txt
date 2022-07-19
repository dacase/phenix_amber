###########################################################################
# The following lines specify the REDFT/RODFT/DHT sizes for which to generate
# specialized codelets.  Currently, only REDFT01/10 of size 8 (used in JPEG).

# e<a><b>_<n> is a hard-coded REDFT<a><b> FFT (DCT) of size <n>
set(E01 e01_8.c)
set(E10 e10_8.c)

add_library(rdft_scalar_r2r OBJECT ${E01} ${E10} codlist.c)
make_pic_if_needed(rdft_scalar_r2r)