CXX=g++ -std=c++11
CFLAGS=-I./include
LDFLAGS=-lopencv_core -lopencv_highgui -lopencv_imgproc

SRC=./src
INCLUDE=./include
BUILD_DIR=./build
BIN_DIR=./build/bin
TOOL_SRC=./tools

MKDIR_P=mkdir -p

all: build_folder $(BUILD_DIR)/RegionGrowth.o $(BUILD_DIR)/TextDetection.o $(BUILD_DIR)/CharRegion.o $(BIN_DIR)/test_main

build_folder:
	$(MKDIR_P) $(BUILD_DIR)
	$(MKDIR_P) $(BIN_DIR)

$(BUILD_DIR)/RegionGrowth.o: $(SRC)/RegionGrowth.cpp
	$(CXX) $^ $(CFLAGS) -c -o $@

$(BUILD_DIR)/TextDetection.o: $(SRC)/TextDetection.cpp
	$(CXX) $^ $(CFLAGS) -c -o $@

$(BUILD_DIR)/CharRegion.o: $(SRC)/CharRegion.cpp
	$(CXX) $^ $(CFLAGS) -c -o $@

$(BIN_DIR)/test_main: $(TOOL_SRC)/test_main.cpp $(BUILD_DIR)/RegionGrowth.o $(BUILD_DIR)/TextDetection.o $(BUILD_DIR)/CharRegion.o
	$(CXX) $^ $(CFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf $(BUILD_DIR)