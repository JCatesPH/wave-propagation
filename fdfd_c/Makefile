# ----------  RULE 1  -------------
# Add source files
CFILES := fdfd_per.c
CDEPS  := 
#OBJFILES := fsource.o
HFILES := domain.h 
CFLAGS := -lm -Wall -lgsl -lgslcblas -std=c99 -I/usr/local/include/gsl -L/usr/local/lib


test : $(CFILES) $(CDEPS) $(HFILES) 
	gcc -o fdfdPerNL_test $(CFILES) $(CFLAGS)

clean:
	rm fdfdPerNL_test