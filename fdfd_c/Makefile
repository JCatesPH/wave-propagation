# ----------  RULE 1  -------------
# Add source files
CFILES := fdfd_per.c
CDEPS  := 
#OBJFILES := fsource.o
HFILES := domain.h
CFLAGS := -lm -Wall -lgsl -lgslcblas -std=c11 -I/usr/local/include/gsl:/home/jmcates/SuiteSparse/include -L/usr/local/lib:/home/jmcates/SuiteSparse/lib


test : $(CFILES) $(CDEPS) $(HFILES) 
	cc -o fdfdPerNL_test $(CFILES) $(CFLAGS)

clean:
	rm fdfdPerNL_test