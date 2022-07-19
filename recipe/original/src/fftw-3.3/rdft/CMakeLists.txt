set(RDFT_SOURCES buffered2.c direct2.c nop2.c rank0-rdft2.c rank-geq2-rdft2.c	
	plan2.c problem2.c solve2.c vrank-geq1-rdft2.c rdft2-rdft.c		
	rdft2-tensor-max-index.c rdft2-inplace-strides.c rdft2-strides.c	
	khc2c.c ct-hc2c.c ct-hc2c-direct.c hc2hc.c dft-r2hc.c dht-r2hc.c dht-rader.c	
	buffered.c conf.c direct-r2r.c direct-r2c.c generic.c	
	hc2hc-direct.c hc2hc-generic.c khc2hc.c kr2c.c kr2r.c indirect.c nop.c	
	plan.c problem.c rank0.c rank-geq2.c rdft-dht.c solve.c		
	vrank-geq1.c vrank3-transpose.c)

add_library(rdft OBJECT ${RDFT_SOURCES})
make_pic_if_needed(rdft)

add_subdirectory(scalar)
add_subdirectory(simd)