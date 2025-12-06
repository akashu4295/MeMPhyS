# --- Compiler Selection ---
# Use: make            -> builds with gcc
#      make GPU=1      -> builds with nvcc

ifeq ($(GPU),1)
    CC = nvcc
    CFLAGS =
    GPU_MSG = "Compiling with NVCC (GPU mode)"
else
    CC = gcc
    CFLAGS = -lm
    GPU_MSG = "Compiling with GCC (CPU mode)"
endif

# Directories
SRC_DIR = header_files
INIT = init/init_TC.c
MAIN = mg_NS_solver.c

# Collect all .c files
SRC = $(wildcard $(SRC_DIR)/*.c)

# Output executable
TARGET = a.out

# Build rule
$(TARGET): $(SRC) $(INIT) $(MAIN)
	@echo $(GPU_MSG)
	$(CC) $(SRC) $(INIT) $(MAIN) $(CFLAGS) -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)
