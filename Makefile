# Compiler and flags
CC = gcc
CFLAGS = -lm

# Directories
SRC_DIR = header_files
INIT = init/init_TC.c
MAIN = mg_NS_solver.c

# Collect all .c files from header_files
SRC = $(wildcard $(SRC_DIR)/*.c)

# Output executable
TARGET = a.out

# Build rule
$(TARGET): $(SRC) $(INIT) $(MAIN)
	$(CC) $(SRC) $(INIT) $(MAIN) $(CFLAGS) -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)
