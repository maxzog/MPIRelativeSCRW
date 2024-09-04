# Compiler
FC = mpif90

# Compiler flags
FFLAGS = -fallow-argument-mismatch -O3 -ffast-math -march=native

# Directories
SRC_DIR = src
OBJ_DIR = obj
CASES_DIR = cases

# Source files
TRACER_SRC = $(SRC_DIR)/tracer_class.f90
TRACER_MAIN = $(CASES_DIR)/tracer.f90

INERTIAL_SRC = $(SRC_DIR)/inertial_class.f90
INERTIAL_MAIN = $(CASES_DIR)/inertial.f90

# Executables
TRACER_EXE = program_tracer
INERTIAL_EXE = program_inertial

directories:
	@mkdir -p ./bin
	@mkdir -p ./obj
	@mkdir -p ./outs

# Default target
all: directories tracer inertial

# Build tracer
tracer: directories $(OBJ_DIR) $(TRACER_SRC) $(TRACER_MAIN)
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -o $(TRACER_EXE) $(TRACER_SRC) $(TRACER_MAIN)

# Build inertial
inertial: directories $(OBJ_DIR) $(INERTIAL_SRC) $(INERTIAL_MAIN)
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -o $(INERTIAL_EXE) $(INERTIAL_SRC) $(INERTIAL_MAIN)

# Clean up
clean:
	rm -f $(TRACER_EXE) $(INERTIAL_EXE)
	rm -r ./obj
