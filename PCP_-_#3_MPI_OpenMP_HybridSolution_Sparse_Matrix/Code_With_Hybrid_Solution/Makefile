################################################################################
# Makefile for PCP 
#
# by Carlos Sá e Ana Sousa 
################################################################################

SHELL = /bin/sh

MPICXX = mpicc
LD  = mpicc

MPI = tp2_hybrid
BIN = bin
BIN_MPI = tp2_hybrid

CXXFLAGS = -O3 -Wall -Wextra -std=c11 -fopenmp
#LIBS = -L/share/apps/mpiP-3.4.1/ -lmpiP -lbfd -liberty -lunwind -lm

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
	$(MPICXX) -M  $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: %.c
	$(MPICXX) -c  $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BIN_DIR)/$(BIN_MPI): $(BUILD_DIR)/$(MPI).o $(BUILD_DIR)/$(MPI).d 
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(MPI).o

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p timing
	@mkdir -p bin

all: checkdirs  $(BIN_DIR)/$(BIN_MPI) 

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/*
