# ----------  RULE 1  -------------
# Add source files
CFILES := fdfd_per_CSparse.c 
CDEPS  := 
#OBJFILES := fsource.o
HFILES := include/domain.h include/cs.h
CFLAGS := -lm -Wall -lrt -lgsl -std=c11 -I./include -L/home/jmcates/SuiteSparse/lib -fexceptions -fPIC -fopenmp -lcxsparse


test : $(CFILES) $(CDEPS) $(HFILES) 
	cc -o fdfdPerNL_test $(CFILES) $(CFLAGS)

clean:
	rm fdfdPerNL_test
