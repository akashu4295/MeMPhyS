############################################
# Compiler Mode Selection
############################################

# ----- CUDA mode (nvcc) -----
ifeq ($(CUDA),1)
    CC = nvcc
    MODE = cuda
    CFLAGS = -O3          # Minimal safe flags for nvcc
    GPU_MSG = "Compiling with NVCC (CUDA mode)"
endif

# ----- OpenACC mode (NVIDIA HPC SDK) -----
ifeq ($(OPENACC),1)
    CC = nvc
    MODE = openacc
    CFLAGS = -acc -Minfo=accel
    GPU_MSG = "Compiling with NVC (OpenACC GPU mode)"
endif

# ----- Default CPU mode (gcc) -----
ifeq ($(MODE),)
    CC = gcc
    MODE = cpu
    GPU_MSG = "Compiling with GCC (CPU mode)"
    CFLAGS = -lm -O3
endif


############################################
# Optional Features (added only if valid)
############################################

# ------- Warnings (only for gcc/nvc, NOT nvcc) ------
ifeq ($(WARN),1)
    ifneq ($(MODE),cuda)
        CFLAGS += -Wall -Wextra -Wshadow -Wpointer-arith
    endif
endif

# ------- Optimisation flags -------
ifeq ($(OPT),1)
    ifeq ($(MODE),cpu)
        CFLAGS += -O3 -march=native
    endif
    ifeq ($(MODE),openacc)
        CFLAGS += -fast
    endif
    # nvcc already has -O3 by default (added above)
endif

# ------- Debugging -------
ifeq ($(DEBUG),1)
    CFLAGS += -g
endif

# ------- Unified Memory (OpenACC or CUDA only) -------
ifeq ($(UM),1)
    ifeq ($(MODE),openacc)
        CFLAGS += -gpu=managed
    endif
    ifeq ($(MODE),cuda)
        CFLAGS += -Xcompiler -rdc=true -arch=sm_80
    endif
endif


############################################
# Directories and Files
############################################

SRC_DIR = header_files
INIT = init/init_TC.c
MAIN = mg_NS_solver.c
SRC = $(wildcard $(SRC_DIR)/*.c)

TARGET = a.out
MODE_FILE = .build_mode


############################################
# Build Rules
############################################

$(TARGET): $(SRC) $(INIT) $(MAIN) $(MODE_FILE)
	@echo $(GPU_MSG)
	$(CC) $(SRC) $(INIT) $(MAIN) $(CFLAGS) -o $(TARGET)

$(MODE_FILE):
	@echo $(MODE) > $(MODE_FILE)

ifeq ($(MODE),$(shell cat $(MODE_FILE) 2>/dev/null))
else
$(TARGET): clean $(MODE_FILE)
endif

clean:
	rm -f $(TARGET) $(MODE_FILE)
