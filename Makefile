# Compiler and flags
FC = mpif90
FLAGS = -fallow-argument-mismatch -O3 -ffast-math -march=native -fopenmp -J$(OBJ_DIR)

# Directory paths
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
OUT_DIR = ./outs
CASES_DIR = ./cases

# Source files
MODULE_SOURCES = $(SRC_DIR)/simulation_class.f90
MODULE_OBJECTS = $(MODULE_SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
SOURCES = $(filter-out $(MODULE_SOURCES), $(wildcard $(SRC_DIR)/*.f90))
OBJECTS = $(SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

# Default target
all: test

# Ensure necessary directories are available
directories:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OUT_DIR)

# Rule to compile .f90 files to .o files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | directories
	$(FC) $(FLAGS) -c $< -o $@

# Test target
test: directories $(BIN_DIR)/program_test

# Link all necessary objects and test case to create the executable
$(BIN_DIR)/program_test: $(CASES_DIR)/test.f90 $(MODULE_OBJECTS) $(OBJECTS)
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

# Clean target
.PHONY: clean
clean:
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(BIN_DIR)/program_test
