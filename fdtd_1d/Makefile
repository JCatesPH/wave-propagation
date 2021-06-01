# --- Variable definitions ---
CC = g++
CFLAGS = -Wall -std=c++17 -I./include 
LDFLAGS = -lm -lrt -fopenmp -fexceptions
OBJFILES = fdtd.o helpers.o
HFILES = include/helpers.h include/parameters.h
TARGET = test

# --- Rule 1 ---
all: clean $(TARGET)

$(TARGET): $(OBJFILES) $(HFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
