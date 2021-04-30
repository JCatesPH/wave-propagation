# ----------  RULE 1  -------------
# Add source files
CFILES := fdtd_2d_PML.c
CDEPS  :=
HFILES := 
CFLAGS := -lm -Wall

pml : $(CFILES) $(CDEPS) $(HFILES)
	g++ -o fdtd2dpml $(CFILES) $(CFLAGS)


# ----------  RULE 2  -------------
# Add source files
CFILES := fdtd_2d_sftf.c
CDEPS  :=
HFILES := 
CFLAGS := -lm -Wall

# Rules and targets
sftf : $(CFILES) $(CDEPS) $(HFILES)
	g++ -o fdtd2dsftf $(CFILES) $(CFLAGS)

clean:
	rm fdtd2dsftf fdtd2dpml