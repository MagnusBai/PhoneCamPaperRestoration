#include <iostream>
#include <opencv2/opencv.hpp>

#include "TextDetection.h"

using namespace std;
using namespace cv;
using namespace DetectText;

Mat loadByteImage(const char * name) {
    Mat image = imread(name);

    if (image.empty()) {
        return Mat();
    }
    cvtColor(image, image, CV_BGR2RGB);
    return image;
}

int mainTextDetection(int argc, char** argv) {
    Mat byteQueryImage = loadByteImage(argv[1]);
    if (byteQueryImage.empty()) {
        cerr << "couldn't load query image" << endl;
        return -1;
    }

    // Detect text in the image
    // Mat output = textDetection(byteQueryImage, atoi(argv[3]));
    Mat edge_swt = swtFilterEdges(byteQueryImage);

    // genThetaHistograms(edge_swt);

    // imwrite(argv[2], output);
    return 0;
}

int main(int argc, char** argv) {
    if ((argc != 2)) {
        cerr << "usage: " << argv[0] << " image_path" << endl;
        return -1;
    }
    return mainTextDetection(argc, argv);
}
