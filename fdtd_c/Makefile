# --- Variable definitions ---
CC = g++
CFLAGS = -Wall -std=c++17 -I./include 
LDFLAGS = -lm -lrt -fopenmp -fexceptions
OBJFILES = fdtd_2D.o helpers.o
HFILES = include/helpers.h include/parameters.h
TARGET = fdtd

# --- Rule 1 ---
all: $(TARGET)

$(TARGET): $(OBJFILES) $(HFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
