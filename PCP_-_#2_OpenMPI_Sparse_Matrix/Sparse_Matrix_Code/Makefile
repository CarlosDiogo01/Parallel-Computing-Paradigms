################################################################################
# Makefile for PCP 
#
# by (CPD-Minho)
################################################################################

SHELL = /bin/sh

CXX = gcc
LD  = gcc

SEQ = tp1_seq
PAR = tp1_paralel
BIN = bin
BIN_SEQ = tp1_seq
BIN_PAR = tp1_paralel

CXXFLAGS   = -O3 -Wall -Wextra -std=c99 -fopenmp

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
	$(CXX) -M $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: %.c
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BIN_DIR)/$(BIN_SEQ): $(BUILD_DIR)/$(SEQ).o $(BUILD_DIR)/$(SEQ).d 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(SEQ).o 

$(BIN_DIR)/$(BIN_PAR): $(BUILD_DIR)/$(PAR).o $(BUILD_DIR)/$(PAR).d 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(PAR).o 

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p timing
	@mkdir -p bin

all: checkdirs  $(BIN_DIR)/$(BIN_SEQ) $(BIN_DIR)/$(BIN_PAR)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/*
