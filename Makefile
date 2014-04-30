all: nm

nm: nulls.c
	gcc nulls.c -o nm -lgsl -lgslcblas -O3 -DHAVE_INLINE
