# --- Variable definitions ---
CC = mpicc
CPPFLAGS = -Wall -std=c++17 -I./include -I/home/jmcates/petsc/include -I/home/jmcates/petsc/arch-linux-c-debug/include -fPIC -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden #-g3
LDFLAGS = -L/home/jmcates/petsc/arch-linux-c-debug/lib -Wl,-rpath,/home/jmcates/petsc/arch-linux-c-debug/lib -lm -lpetsc -lstdc++
OBJFILES = util_petsc.o fdfd_per_petsc.o
TARGET = test

# --- Rule 1 ---
all: $(TARGET)

$(TARGET): $(OBJFILES) include/domain.h
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
