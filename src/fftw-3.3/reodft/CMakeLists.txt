set(REODFT_SOURCES conf.c reodft.h reodft010e-r2hc.c	
	reodft11e-radix2.c reodft11e-r2hc-odd.c redft00e-r2hc-pad.c	
	rodft00e-r2hc-pad.c reodft00e-splitradix.c)

add_library(reodft OBJECT ${REODFT_SOURCES})
make_pic_if_needed(reodft)