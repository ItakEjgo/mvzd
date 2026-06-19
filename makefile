CXX := /usr/bin/g++-14

CXXFLAGS := -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++20 -Wall -Iinclude
LDFLAGS := -L/usr/local/lib -Wl,-R/usr/local/lib
LIBS := -ljemalloc

SRC_DIR := src
BUILD_DIR := build

all: $(BUILD_DIR)/main

$(BUILD_DIR)/main: $(SRC_DIR)/main.cpp $(SRC_DIR)/cpamz.hpp $(SRC_DIR)/mvq.hpp $(SRC_DIR)/hilbert.c $(SRC_DIR)/hilbert.h $(SRC_DIR)/cpambb.hpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/main $(SRC_DIR)/main.cpp $(SRC_DIR)/hilbert.c $(LDFLAGS) $(LIBS)

clean:
	rm -rf $(BUILD_DIR)
