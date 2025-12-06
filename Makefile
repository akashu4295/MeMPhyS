# --- Compiler Selection ---
ifeq ($(GPU),1)
    CC = nvcc
    CFLAGS =
    MODE = gpu
    GPU_MSG = "Compiling with NVCC (GPU mode)"
else
    CC = gcc
    CFLAGS = -lm
    MODE = cpu
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
MODE_FILE = .build_mode

# Build rule
$(TARGET): $(SRC) $(INIT) $(MAIN) $(MODE_FILE)
	@echo $(GPU_MSG)
	$(CC) $(SRC) $(INIT) $(MAIN) $(CFLAGS) -o $(TARGET)

# Track build mode so switching triggers rebuild
$(MODE_FILE):
	@echo $(MODE) > $(MODE_FILE)

# Rebuild if mode changes
ifeq ($(MODE),$(shell cat $(MODE_FILE) 2>/dev/null))
else
$(TARGET): clean $(MODE_FILE)
endif

# Clean rule
clean:
	rm -f $(TARGET) $(MODE_FILE)
