CXX=g++ -std=c++11

LDFLAGS=-static-libstdc++
CFLAGS=-I./include
OPENCV_INCL=
OPENCV_LIBS=-lopencv_core -lopencv_highgui -lopencv_imgproc
MRPT_INCL=-I/usr/include/eigen3 -I/usr/include/mrpt/mrpt-config -I/usr/include/suitesparse -I/usr/include/mrpt/bayes/include -I/usr/include/mrpt/graphs/include -I/usr/include/mrpt/vision/include -I/usr/include/mrpt/tfest/include -I/usr/include/mrpt/maps/include -I/usr/include/mrpt/obs/include -I/usr/include/mrpt/opengl/include -I/usr/include/mrpt/base/include -I/usr/include/mrpt/slam/include -I/usr/include/mrpt/gui/include -I/usr/include/mrpt/topography/include
MRPT_LIBS=-lmrpt-obs -lmrpt-slam -lmrpt-opengl -lmrpt-base -lmrpt-gui

SRC=./src
INCLUDE=./include
BUILD_DIR=./build
BIN_DIR=./build/bin
TOOL_SRC=./tools

MKDIR_P=mkdir -p

all: build_folder $(BIN_DIR)/test_main

build_folder:
	$(MKDIR_P) $(BUILD_DIR)
	$(MKDIR_P) $(BIN_DIR)

$(BUILD_DIR)/RegionGrowth.o: $(SRC)/RegionGrowth.cpp
	$(CXX) $^ $(CFLAGS) $(OPENCV_INCL) -c -o $@

$(BUILD_DIR)/TextDetection.o: $(SRC)/TextDetection.cpp
	$(CXX) $^ $(CFLAGS) $(OPENCV_INCL) $(MRPT_INCL) -c -o $@

$(BUILD_DIR)/CharRegion.o: $(SRC)/CharRegion.cpp
	$(CXX) $^ $(CFLAGS) $(OPENCV_INCL) $(MRPT_INCL) -c -o $@

# $(BIN_DIR)/test_main: $(TOOL_SRC)/test_main.cpp $(BUILD_DIR)/RegionGrowth.o $(BUILD_DIR)/TextDetection.o $(BUILD_DIR)/CharRegion.o
# 	$(CXX) $^ $(LDFLAGS) $(CFLAGS) $(OPENCV_INCL) $(OPENCV_LIBS) $(MRPT_LIBS) -o $@

$(BUILD_DIR)/test_main.o: $(TOOL_SRC)/test_main.cpp
	$(CXX) $^ $(CFLAGS) $(OPENCV_INCL) -c -o $@

$(BUILD_DIR)/PagePartition.o : $(SRC)/PagePartition.cpp
	$(CXX) $^ $(CFLAGS) -c -o $@

$(BIN_DIR)/test_main: $(BUILD_DIR)/test_main.o $(BUILD_DIR)/CharRegion.o $(BUILD_DIR)/RegionGrowth.o $(BUILD_DIR)/TextDetection.o $(BUILD_DIR)/PagePartition.o
	$(CXX) $^ $(LDFLAGS) $(CFLAGS) $(OPENCV_INCL) $(OPENCV_LIBS) $(MRPT_LIBS) -o $@

clean:
	rm -rf $(BUILD_DIR)

clean-im:
	rm *.png
