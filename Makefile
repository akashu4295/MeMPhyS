# Compiler Mode Selection
# Usage:
#   make OPENACC=1    # for OpenACC GPU mode
#   make DEBUG=1      # to include debug symbols
#   make              # for default CPU mode
# ----- OpenACC mode -----
ifeq ($(OPENACC),1)
    CC = nvc
    MODE = openacc
    CFLAGS = -acc -Minfo=accel -O3 -Wall -Wextra -Wshadow -Wpointer-arith
    GPU_MSG = "Compiling with NVC (OpenACC GPU mode)"
endif

# ----- Default CPU mode -----
ifeq ($(MODE),)
    CC = gcc
    MODE = cpu
    CFLAGS = -O3 -march=native -Wall -Wextra -Wshadow -Wpointer-arith
    LDFLAGS = -lm
    GPU_MSG = "Compiling with GCC (CPU mode)"
endif


# Setting Debug flag

ifeq ($(DEBUG),1)
    CFLAGS += -g
endif


# Files

SRC_DIR = header_files
SRC = $(wildcard $(SRC_DIR)/*.c)
INIT = init.c
MAIN = mg_NS_solver.c
TARGET = a.out
MODE_FILE = .build_mode


# Build Rules

all: check_mode $(TARGET)

$(TARGET): $(SRC) $(INIT) $(MAIN)
	@echo $(GPU_MSG)
	$(CC) $(SRC) $(INIT) $(MAIN) $(CFLAGS) -o $(TARGET) $(LDFLAGS)

check_mode:
	@if [ -f $(MODE_FILE) ] && [ "`cat $(MODE_FILE)`" != "$(MODE)" ]; then \
	    echo "Build mode changed → full rebuild"; \
	    rm -f $(TARGET); \
	fi
	@echo $(MODE) > $(MODE_FILE)

clean:
	rm -f $(TARGET) $(MODE_FILE)