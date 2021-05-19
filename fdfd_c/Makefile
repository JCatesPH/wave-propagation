# --- Variable definitions ---
CC = gcc
CFLAGS = -Wall -std=c11 -I./include -L/home/jmcates/SuiteSparse/lib 
LDFLAGS = -lm -lrt -lgsl -fopenmp -lcxsparse -fexceptions -fPIC 
OBJFILES = util.o fdfd_per_CSparse.o
TARGET = fdfd_per

# --- Rule 1 ---
all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
