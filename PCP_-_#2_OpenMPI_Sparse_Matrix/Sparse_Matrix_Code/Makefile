################################################################################
# Makefile for PCP 
#
# by Carlos Sá e Ana Sousa 
################################################################################

SHELL = /bin/sh

CXX = mpicc
LD  = mpicc

MPI = tp2_mpi
BIN = bin
BIN_MPI = tp2_mpi

SRC_DIR = src
BIN_DIR = bin
BUILD_DIR = build
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst src/*.c,build/*.o,$(SRC))
DEPS = $(patsubst build/*.o,build/*.d,$(OBJ))

vpath %.c $(SRC_DIR)


################################################################################
# Rules
################################################################################

.DEFAULT_GOAL = all


$(BUILD_DIR)/%.d: %.c
	$(CXX) -M  $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: %.c
	$(CXX) -c  $(INCLUDES) $< -o $@

$(BIN_DIR)/$(BIN_MPI): $(BUILD_DIR)/$(MPI).o $(BUILD_DIR)/$(MPI).d 
	$(CXX) -o $@ $(BUILD_DIR)/$(MPI).o 

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p bin

all: checkdirs  $(BIN_DIR)/$(BIN_MPI) 

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/*
