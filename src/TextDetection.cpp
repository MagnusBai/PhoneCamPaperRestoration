/*
 Copyright 2012 Andrew Perrault and Saurav Kumar.

 This file is part of DetectText.

 DetectText is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 DetectText is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with DetectText.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/unordered_map.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
/*#include <graph/adjacency_list.hpp>
 #include <graph/graph_traits.hpp>
 #include <graph/connected_components.hpp>
 #include <unordered_map.hpp>
 #include <graph/floyd_warshall_shortest.hpp>
 #include <numeric/ublas/matrix.hpp>
 #include <numeric/ublas/io.hpp> */
#include <cassert>
#include <cmath>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include <time.h>
#include <utility>
#include <algorithm>
#include <vector>

#include <limits>

//#include <opencv/cv.h>
//#include <opencv/cxcore.h>
//#include <opencv/highgui.h>
#include <unordered_map>
#include <string>



#include "TextDetection.h"

#include "CharRegion.h"

// #include <mrpt/math/lightweight_geom_data.h>		// mrpt::math::TPoint2D TLine2D
#include <mrpt/math/ransac_applications.h>
#include <mrpt/gui/CDisplayWindow3D.h>
#include <mrpt/gui/CDisplayWindowPlots.h>
#include <mrpt/random.h>
#include <mrpt/utils/CTicTac.h>
#include <mrpt/utils/metaprogramming.h>
#include <mrpt/poses/CPose3D.h>
#include <mrpt/opengl/CGridPlaneXY.h>
#include <mrpt/opengl/CPointCloud.h>
#include <mrpt/opengl/stock_objects.h>
#include <mrpt/opengl/CTexturedPlane.h>

#include "PagePartition.h"

using namespace cv;
using namespace std;

#define PI CV_PI

namespace DetectText {

const Scalar BLUE(255, 0, 0);
const Scalar GREEN(0, 255, 0);
const Scalar RED(0, 0, 255);

template<typename T>
void clac1DHistogram(const vector<T>& data, const int count_hist,
		const pair<T, T>& range, vector<int>& hist,
		vector<vector<int>>& hist_ids) {
	hist = vector<int>(count_hist, 0);
	hist_ids = vector<vector<int>>(count_hist, vector<int>(0));
	T cell_range = (range.second + static_cast<T>(0.001f) - range.first)
			/ static_cast<T>(count_hist);
	for (int i = 0; i < data.size(); ++i) {
		int id = static_cast<int>(floor(data[i] / cell_range));
		++hist[id];
		hist_ids[id].push_back(i);
	}
}

template<typename T>
void conv1DSignal(const vector<T>& signal, const vector<float>& convFilter,
		vector<T>& res_signal) {
	res_signal = vector<T>(signal.size());	// preparation for result distribution

	int conv_size = convFilter.size();
	assert(conv_size % 2 == 1);	// conv filter's receptive field should be odd
	int pos_mid = conv_size / 2;
	// take 7 sized conv filter as example
	float mid_multiplier = convFilter[pos_mid];				// is at filter[3]
	vector<float> prefix_multipliers(pos_mid);	// filter[2: -1: -1]
	vector<float> suffix_multipliers(pos_mid);	// filter[4: 7]
	for (int i = 0; i < pos_mid; ++i) {
		prefix_multipliers[i] = convFilter[pos_mid - 1 - i];
		suffix_multipliers[i] = convFilter[pos_mid + 1 + i];
	}

	int i_pos = 0;
	for (typename vector<T>::const_iterator it_mid = signal.begin(); it_mid < signal.end();
			++it_mid) {
		float conv_res = static_cast<float>(*it_mid) * mid_multiplier;
		// go through prefix
		for (int i = 0; i < pos_mid; ++i) {
			typename vector<T>::const_iterator it = it_mid;
			if (it == signal.begin()) {
				it = signal.end();
			}
			it--;
			conv_res += static_cast<float>(*it) * prefix_multipliers[i];
		}
		// go through suffix
		for (int i = 0; i < pos_mid; ++i) {
			typename vector<T>::const_iterator it = it_mid;
			it++;
			if (it == signal.end()) {
				it = signal.begin();
			}
			conv_res += static_cast<float>(*it) * suffix_multipliers[i];
		}

		res_signal[i_pos] = static_cast<T>(round(conv_res));
		++i_pos;
	}
}

std::vector<SWTPointPair2i> findBoundingBoxes(
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<Chain> & chains, std::vector<SWTPointPair2d> & compBB,
		Mat& output) {
	std::vector<SWTPointPair2i> bb;
	bb.reserve(chains.size());
	for (auto& chainit : chains) {
		int minx = output.cols;
		int miny = output.rows;
		int maxx = 0;
		int maxy = 0;
		for (std::vector<int>::const_iterator cit = chainit.components.begin();
				cit != chainit.components.end(); cit++) {
			miny = std::min(miny, compBB[*cit].first.y);
			minx = std::min(minx, compBB[*cit].first.x);
			maxy = std::max(maxy, compBB[*cit].second.y);
			maxx = std::max(maxx, compBB[*cit].second.x);
		}
		Point2i p0(minx, miny);
		Point2i p1(maxx, maxy);
		SWTPointPair2i pair(p0, p1);
		bb.push_back(pair);
	}
	return bb;
}

std::vector<SWTPointPair2i> findBoundingBoxes(
		std::vector<std::vector<SWTPoint2d> > & components, Mat& output) {
	std::vector<SWTPointPair2i> bb;
	bb.reserve(components.size());
	for (auto& compit : components) {
		int minx = output.cols;
		int miny = output.rows;
		int maxx = 0;
		int maxy = 0;
		for (auto& it : compit) {
			miny = std::min(miny, it.y);
			minx = std::min(minx, it.x);
			maxy = std::max(maxy, it.y);
			maxx = std::max(maxx, it.x);
		}
		Point2i p0(minx, miny);
		Point2i p1(maxx, maxy);
		SWTPointPair2i pair(p0, p1);
		bb.push_back(pair);
	}
	return bb;
}

void normalizeImage(const Mat& input, Mat& output) {
	assert(input.depth() == CV_32F);
	assert(input.channels() == 1);
	assert(output.depth() == CV_32F);
	assert(output.channels() == 1);

	float maxVal = 0;
	float minVal = 1e100;
	for (int row = 0; row < input.rows; row++) {
		const float* ptr = (const float*) input.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			if (*ptr < 0) {
			} else {
				maxVal = std::max(*ptr, maxVal);
				minVal = std::min(*ptr, minVal);
			}
			ptr++;
		}
	}

	float difference = maxVal - minVal;
	for (int row = 0; row < input.rows; row++) {
		const float* ptrin = (const float*) input.ptr(row);
		float* ptrout = (float*) output.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			if (*ptrin < 0) {
				*ptrout = 1;
			} else {
				*ptrout = ((*ptrin) - minVal) / difference;
			}
			ptrout++;
			ptrin++;
		}
	}
}

void renderComponents(const Mat& SWTImage,
		std::vector<std::vector<SWTPoint2d> > & components, Mat& output) {
	output.setTo(0);

	for (auto& component : components) {
		for (auto& pit : component) {
			output.at<float>(pit.y, pit.x) = SWTImage.at<float>(pit.y, pit.x);
		}
	}
	for (int row = 0; row < output.rows; row++) {
		float* ptr = (float*) output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr == 0) {
				*ptr = -1;
			}
			ptr++;
		}
	}
	float maxVal = 0;
	float minVal = 1e100;
	for (int row = 0; row < output.rows; row++) {
		const float* ptr = (const float*) output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr == 0) {
			} else {
				maxVal = std::max(*ptr, maxVal);
				minVal = std::min(*ptr, minVal);
			}
			ptr++;
		}
	}
	float difference = maxVal - minVal;
	for (int row = 0; row < output.rows; row++) {
		float* ptr = (float*) output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr < 1) {
				*ptr = 1;
			} else {
				*ptr = ((*ptr) - minVal) / difference;
			}
			ptr++;
		}
	}

}

void renderComponentsWithBoxes(Mat& SWTImage,
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<SWTPointPair2d> & compBB, Mat& output) {
	Mat outTemp(output.size(), CV_32FC1);

	renderComponents(SWTImage, components, outTemp);

	std::vector<SWTPointPair2i> bb;
	bb.reserve(compBB.size());
	for (auto& it : compBB) {
		Point2i p0 = cvPoint(it.first.x, it.first.y);
		Point2i p1 = cvPoint(it.second.x, it.second.y);
		SWTPointPair2i pair(p0, p1);
		bb.push_back(pair);
	}

	Mat out(output.size(), CV_8UC1);
	outTemp.convertTo(out, CV_8UC1, 255.);
	cvtColor(out, output, CV_GRAY2RGB);

	int count = 0;
	for (auto it : bb) {
		Scalar c;
		if (count % 3 == 0) {
			c = BLUE;
		} else if (count % 3 == 1) {
			c = GREEN;
		} else {
			c = RED;
		}
		count++;
		rectangle(output, it.first, it.second, c, 2);
	}
}

void renderChainsWithBoxes(Mat& SWTImage,
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<Chain> & chains, std::vector<SWTPointPair2d> & compBB,
		Mat& output) {
	// keep track of included components
	std::vector<bool> included;
	included.reserve(components.size());
	for (unsigned int i = 0; i != components.size(); i++) {
		included.push_back(false);
	}
	for (Chain& it : chains) {
		for (std::vector<int>::iterator cit = it.components.begin();
				cit != it.components.end(); cit++) {
			included[*cit] = true;
		}
	}
	std::vector < std::vector<SWTPoint2d> > componentsRed;
	for (unsigned int i = 0; i != components.size(); i++) {
		if (included[i]) {
			componentsRed.push_back(components[i]);
		}
	}
	Mat outTemp(output.size(), CV_32FC1);

	std::cout << componentsRed.size() << " components after chaining"
			<< std::endl;
	renderComponents(SWTImage, componentsRed, outTemp);
	std::vector<SWTPointPair2i> bb;
	bb = findBoundingBoxes(components, chains, compBB, outTemp);

	Mat out(output.size(), CV_8UC1);
	outTemp.convertTo(out, CV_8UC1, 255);
	cvtColor(out, output, CV_GRAY2RGB);

	int count = 0;
	for (auto& it : bb) {
		CvScalar c;
		if (count % 3 == 0) {
			c = BLUE;
		} else if (count % 3 == 1) {
			c = GREEN;
		} else {
			c = RED;
		}
		count++;
		rectangle(output, it.first, it.second, c, 2);
	}
}

void renderChains(Mat& SWTImage,
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<Chain> & chains, Mat& output) {
	// keep track of included components
	std::vector<bool> included;
	included.reserve(components.size());
	for (unsigned int i = 0; i != components.size(); i++) {
		included.push_back(false);
	}
	for (std::vector<Chain>::iterator it = chains.begin(); it != chains.end();
			it++) {
		for (std::vector<int>::iterator cit = it->components.begin();
				cit != it->components.end(); cit++) {
			included[*cit] = true;
		}
	}
	std::vector < std::vector<SWTPoint2d> > componentsRed;
	for (unsigned int i = 0; i != components.size(); i++) {
		if (included[i]) {
			componentsRed.push_back(components[i]);
		}
	}
	std::cout << componentsRed.size() << " components after chaining"
			<< std::endl;
	Mat outTemp(output.size(), CV_32FC1);
	renderComponents(SWTImage, componentsRed, outTemp);
	outTemp.convertTo(output, CV_8UC1, 255);

}

Mat textDetection(const Mat& input, bool dark_on_light) {
	assert(input.depth() == CV_8U);
	assert(input.channels() == 3);

	std::cout << "Running textDetection with dark_on_light " << dark_on_light
			<< std::endl;

	// Convert to grayscale
	Mat grayImage(input.size(), CV_8UC1);
	cvtColor(input, grayImage, CV_RGB2GRAY);
	// Create Canny Image
	double threshold_low = 175;
	double threshold_high = 320;
	Mat edgeImage(input.size(), CV_8UC1);
	Canny(grayImage, edgeImage, threshold_low, threshold_high, 3);
	imwrite("canny.png", edgeImage);

	// Create gradient X, gradient Y
	Mat gaussianImage(input.size(), CV_32FC1);
	grayImage.convertTo(gaussianImage, CV_32FC1, 1. / 255.);
	GaussianBlur(gaussianImage, gaussianImage, Size(5, 5), 0);
	Mat gradientX(input.size(), CV_32FC1);
	Mat gradientY(input.size(), CV_32FC1);
	Scharr(gaussianImage, gradientX, -1, 1, 0);
	Scharr(gaussianImage, gradientY, -1, 0, 1);
	GaussianBlur(gradientX, gradientX, Size(3, 3), 0);
	GaussianBlur(gradientY, gradientY, Size(3, 3), 0);

	// Calculate SWT and return ray vectors
	std::vector<Ray> rays;
	Mat SWTImage(input.size(), CV_32FC1);
	for (int row = 0; row < input.rows; row++) {
		float* ptr = (float*) SWTImage.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			*ptr++ = -1;
		}
	}
	strokeWidthTransform(edgeImage, gradientX, gradientY, dark_on_light,
			SWTImage, rays);
	SWTMedianFilter(SWTImage, rays);

	Mat output2(input.size(), CV_32FC1);
	normalizeImage(SWTImage, output2);
	Mat saveSWT(input.size(), CV_8UC1);
	output2.convertTo(saveSWT, CV_8UC1, 255);
	imwrite("SWT.png", saveSWT);

	// Calculate legally connect components from SWT and gradient image.
	// return type is a vector of vectors, where each outer vector is a component and
	// the inner vector contains the (y,x) of each pixel in that component.
	std::vector < std::vector<SWTPoint2d> > components =
			findLegallyConnectedComponents(SWTImage, rays);

	// Filter the components
	std::vector < std::vector<SWTPoint2d> > validComponents;
	std::vector<SWTPointPair2d> compBB;
	std::vector<Point2dFloat> compCenters;
	std::vector<float> compMedians;
	std::vector<SWTPoint2d> compDimensions;
	filterComponents(SWTImage, components, validComponents, compCenters,
			compMedians, compDimensions, compBB);

	Mat output3(input.size(), CV_8UC3);
	renderComponentsWithBoxes(SWTImage, validComponents, compBB, output3);
	imwrite("components.png", output3);
	//

	// Make chains of components
	std::vector<Chain> chains;
	chains = makeChains(input, validComponents, compCenters, compMedians,
			compDimensions, compBB);

	Mat output4(input.size(), CV_8UC1);
	renderChains(SWTImage, validComponents, chains, output4);
	imwrite("text.png", output4);

	Mat output5(input.size(), CV_8UC3);
	cvtColor(output4, output5, CV_GRAY2RGB);

	/*IplImage * output =
	 cvCreateImage ( input.size(), CV_8UC3 );
	 renderChainsWithBoxes ( SWTImage, validComponents, chains, compBB, output); */
	return output5;
}

Mat swtFilterEdges(const Mat& input) {
	bool dark_on_light = true;
	assert(input.depth() == CV_8U);
	assert(input.channels() == 3);

	imwrite("original.png", input);

	// Convert to grayscale
	Mat grayImage(input.size(), CV_8UC1);
	cvtColor(input, grayImage, CV_RGB2GRAY);
	// Create Canny Image
	double threshold_low = 175;
	double threshold_high = 320;
	Mat edgeImage(input.size(), CV_8UC1);
	Canny(grayImage, edgeImage, threshold_low, threshold_high, 3);
	imwrite("canny.png", edgeImage);
	// 		edge is repr as 255 in grayImage && bg is repr as 0

	// Create gradient X, gradient Y
	Mat gaussianImage(input.size(), CV_32FC1);
	grayImage.convertTo(gaussianImage, CV_32FC1, 1. / 255.);
	GaussianBlur(gaussianImage, gaussianImage, Size(5, 5), 0);
	Mat gradientX(input.size(), CV_32FC1);
	Mat gradientY(input.size(), CV_32FC1);
	Scharr(gaussianImage, gradientX, -1, 1, 0);
	Scharr(gaussianImage, gradientY, -1, 0, 1);
	GaussianBlur(gradientX, gradientX, Size(3, 3), 0);
	GaussianBlur(gradientY, gradientY, Size(3, 3), 0);

	// Calculate SWT and return ray vectors
	std::vector<Ray> rays;
	Mat SWTImage(input.size(), CV_32FC1);
	for (int row = 0; row < input.rows; row++) {
		float* ptr = (float*) SWTImage.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			*ptr++ = -1;
		}
	}
	strokeWidthTransform(edgeImage, gradientX, gradientY, dark_on_light,
			SWTImage, rays);
	SWTMedianFilter(SWTImage, rays);

//	Mat output2(input.size(), CV_32FC1);
//	normalizeImage(SWTImage, output2);
//	Mat saveSWT(input.size(), CV_8UC1);
//	output2.convertTo(saveSWT, CV_8UC1, 255);
//	imwrite("SWT.png", saveSWT);

	// Calculate legally connect components from SWT and gradient image.
	// return type is a vector of vectors, where each outer vector is a component and
	// the inner vector contains the (y,x) of each pixel in that component.
	std::vector < std::vector<SWTPoint2d> > components =
			findLegallyConnectedComponents(SWTImage, rays);

	// Filter the components
	std::vector < std::vector<SWTPoint2d> > validComponents;
	std::vector<SWTPointPair2d> compBB;
	std::vector<Point2dFloat> compCenters;
	std::vector<float> compMedians;
	std::vector<SWTPoint2d> compDimensions;
	filterComponents(SWTImage, components, validComponents, compCenters,
			compMedians, compDimensions, compBB);

	Mat char_plot(input.size(), CV_8UC1, cv::Scalar(0));
	int plot_scope = 4;
	for (int i = 0; i < validComponents.size(); ++i) {
		vector<SWTPoint2d>& char_vec = validComponents.at(i);
		for (auto& pt : char_vec) {
			int h0_ = pt.y;
			int w0_ = pt.x;
			for (int h_ = h0_ - plot_scope; h_ < h0_ + plot_scope; ++h_) {
				for (int w_ = w0_ - plot_scope; w_ < w0_ + plot_scope; ++w_) {
					if (h_ > 0 && w_ > 0 && h_ < char_plot.rows
							&& w_ < char_plot.cols) {
						char_plot.at < uchar > (h_, w_) = 255;
					}
				}
			}
		}
	}
	imwrite("chars.png", char_plot);

	// SEGMENT EACH PART OF CANDIDATE LETTERS PATCH
	IplImage ipl_char_plot = char_plot;		// maybe import memory leaking...

	//// alloc pres_mat
	int * pres_mat;
	pres_mat = new int[ipl_char_plot.width * ipl_char_plot.height];
	for (int i = 0; i < ipl_char_plot.width * ipl_char_plot.height; i++) {
		pres_mat[i] = INIT;
	}

	//    find init seed point
	int seedx, seedy;
	bool STOP_FLAG = false;
	for (int h = 0; h < char_plot.rows; ++h) {
		for (int w = 0; w < char_plot.cols; ++w) {
			if (char_plot.at < uchar > (h, w) == 0) {
				seedx = w;
				seedy = h;
				STOP_FLAG = true;
				break;
			}
		}
		if (STOP_FLAG) {
			break;
		}
	}

	int count_ids = regionGrowth_pro(&ipl_char_plot, pres_mat,
			ipl_char_plot.width, ipl_char_plot.height, seedx, seedy, 50);

	// CANNY-EDGES in the SWT-WORD MASK
	Mat edge_swt(input.size(), CV_8UC1);
	bitwise_and(char_plot, edgeImage, edge_swt);
	imwrite("edges_swt.png", edge_swt);

	// chars-region vars
	vector<CharRegion> charRegionArray;
	// auxiliary vars
	int ori_width = ipl_char_plot.width;
	int ori_height = ipl_char_plot.height;

	int count_char_regions = count_ids / 2;
	vector<int> x_min(count_char_regions, numeric_limits<int>::max());
	vector<int> y_min(count_char_regions, numeric_limits<int>::max());
	vector<int> x_max(count_char_regions, -1);
	vector<int> y_max(count_char_regions, -1);
	vector<int> count_region_pixs(count_char_regions, 0);
	vector<int> count_region_contours(count_char_regions, 0);
	vector<int> count_char_edge(count_char_regions, 0);
	for (int h = 0; h < ori_height; ++h) {
		for (int w = 0; w < ori_width; ++w) {
			int id0 = pres_mat[h * ori_width + w];
			if (id0 < 1) {	// is not independent region, but is a background
				continue;
			}
			int region_id = (id0 + 1) / 2 - 1;
			bool is_cont = ((id0 % 2) == 0);

			if (w < x_min[region_id]) {
				x_min[region_id] = w;
			}
			if (w > x_max[region_id]) {
				x_max[region_id] = w;
			}
			if (h < y_min[region_id]) {
				y_min[region_id] = h;
			}
			if (h > y_max[region_id]) {
				y_max[region_id] = h;
			}

			count_region_pixs[region_id] += 1;
			if (is_cont) {
				count_region_contours[region_id] += 1;
			}
			if (edge_swt.at < uchar > (h, w) != 0) {
				count_char_edge[region_id] += 1;
			}

		}
	}

	// init charRegionArray
	for (int i = 0; i < count_char_regions; ++i) {
		charRegionArray.push_back(
				CharRegion(count_region_pixs[i], count_region_contours[i],
						count_char_edge[i], x_max[i], x_min[i], y_max[i],
						y_min[i]));
	}

	// assign charRegionArray
	for (int h = 0; h < ori_height; ++h) {
		for (int w = 0; w < ori_width; ++w) {
			int id0 = pres_mat[h * ori_width + w];
			if (id0 < 1) {	// is not independent region, but is a background
				continue;
			}
			int region_id = (id0 + 1) / 2 - 1;
			bool is_cont = ((id0 % 2) == 0);

			charRegionArray[region_id].push_region_pix(w, h);
			if (is_cont) {
				charRegionArray[region_id].push_region_contours(w, h);
			}
			if (edge_swt.at < uchar > (h, w) != 0) {
				charRegionArray[region_id].push_edges(w, h);
			}
		}
	}

	// temp
	string char_filename = "char_img_";
	for (int i = 0; i < count_char_regions; ++i) {
		string name = char_filename + to_string(i) + ".png";
		cout << name << endl;
		if (charRegionArray[i].is_measurable()) {
			cout << "\n\n finding " << i << "lines params \n\n\n";
			charRegionArray[i].ransac_find_lines();
			// charRegionArray[i].plot_char_region(name);
		}
	}

	// TEMP: plot lines
	Mat lines_im(input.size(), CV_8UC3, CharRegion::solarized_palette["base3"]);
	for (int i = 0; i < count_char_regions; ++i) {
		// draw region_pix
		// cout << "------" << charRegionArray[i].count_pixs << endl;

		for (int j = 0; j < charRegionArray[i].count_pixs; ++j) {
			int x = charRegionArray[i].region_pixs[j * 2];
			int y = charRegionArray[i].region_pixs[j * 2 + 1];
			// im.at < uchar > (y, x) = 50;
			uchar bgr[3];
			if (charRegionArray[i].const_nlines > 0) {
				bgr[0] = CharRegion::solarized_palette["yellow"][0];
				bgr[1] = CharRegion::solarized_palette["yellow"][1];
				bgr[2] = CharRegion::solarized_palette["yellow"][2];
			} else {
				bgr[0] = CharRegion::solarized_palette["green"][0];
				bgr[1] = CharRegion::solarized_palette["green"][1];
				bgr[2] = CharRegion::solarized_palette["green"][2];
			}
			lines_im.at < cv::Vec3b > (y, x)[0] = bgr[0];
			lines_im.at < cv::Vec3b > (y, x)[1] = bgr[1];
			lines_im.at < cv::Vec3b > (y, x)[2] = bgr[2];
		}
//		// draw edge
		for (int j = 0; j < charRegionArray[i].count_edges; ++j) {
			int x = charRegionArray[i].char_edge[j * 2];
			int y = charRegionArray[i].char_edge[j * 2 + 1];
			lines_im.at < cv::Vec3b > (y, x)[0] =
					CharRegion::solarized_palette["base03"][0];
			lines_im.at < cv::Vec3b > (y, x)[1] =
					CharRegion::solarized_palette["base03"][1];
			lines_im.at < cv::Vec3b > (y, x)[2] =
					CharRegion::solarized_palette["base03"][2];
		}
		// draw line
		// cout << charRegionArray[i].const_nlines << " : ";
		// cout << charRegionArray[i].line_paramsPts_global.size() << endl ;
		for (int j = 0; j < charRegionArray[i].const_nlines; ++j) {
			if (j >= 3) {
				// break;
			}
			cv::line(lines_im,
					cv::Point(
							charRegionArray[i].line_paramsPts_global[j * 4 + 0],
							charRegionArray[i].line_paramsPts_global[j * 4 + 1]),
					cv::Point(
							charRegionArray[i].line_paramsPts_global[j * 4 + 2],
							charRegionArray[i].line_paramsPts_global[j * 4 + 3]),
					CharRegion::solarized_palette["yellow"], 1);
		}

	}
	imwrite("im_lines.png", lines_im);
	// TEMP: plot lines

	vector<pair<float, float>> angles_scores_set;
	for (int i = 0; i < count_char_regions; ++i) {
		vector<pair<float, float>>& newbies =
				charRegionArray[i].region_angle_scores;
		angles_scores_set.insert(angles_scores_set.end(), newbies.begin(),
				newbies.end());
	}
	cout << endl << endl << endl;
	for (int i = 0; i < angles_scores_set.size(); ++i) {
		cout << "(" << angles_scores_set[i].first << ", "
				<< angles_scores_set[i].second;
		cout << ")\t";
	}
	cout << endl << endl << endl;

	// TEMP: plot region lines
	for (int i = 0; i < count_char_regions; ++i) {
		CharRegion& region = charRegionArray[i];
		if (region.const_nlines > 0) {
			Point pt01 = Point(region.x_min, region.y_min);
			Point pt02 = Point(region.x_max, region.y_max);

			rectangle(lines_im, Rect(pt01, pt02),
					CharRegion::solarized_palette["cyan"], 2);
			for (int j = 0; j < region.line_clusters_ids.size(); ++j) {
				int line_id = region.line_clusters_ids[j][0];
				int xi1, yi1, xi2, yi2;
				xi1 = region.bb_intersection_pts[line_id * 4 + 0];
				yi1 = region.bb_intersection_pts[line_id * 4 + 1];
				xi2 = region.bb_intersection_pts[line_id * 4 + 2];
				yi2 = region.bb_intersection_pts[line_id * 4 + 3];

				if (xi1 != -1) {
					line(lines_im, Point(xi1, yi1), Point(xi2, yi2),
							CharRegion::solarized_palette["cyan"], 1);
				}
			}
		}
	}
	imwrite("region_lines.png", lines_im);
	// TEMP: plot region lines

//	// GET GLOBAL AVERAGE SLOPE
//	// get to wide ranged lines ---- abort this method
//	int count_rslopes = 0;		// rslope is short for region-slope
//	for(int i = 0; i < count_char_regions; ++i) {
//		CharRegion& region = charRegionArray[i];
//		count_rslopes += region.line_clusters_ids.size();
//	}
//	vector<int> region_ids = vector<int>(count_rslopes, 0);
//	vector<float> rslope_angles = vector<float>(count_rslopes, 0.f);
//	vector<float> rslope_scores = vector<float>(count_rslopes, 0.f);
//	int i_rslope = 0;
//	for(int i = 0; i < count_char_regions; ++i) {
//		CharRegion& region = charRegionArray[i];
//		for(int j=0; j<region.region_angle_scores.size(); ++j) {
//			region_ids[i_rslope] = i;
//			rslope_angles[i_rslope] = region.region_angle_scores[j].first;
//			rslope_scores[i_rslope] = region.region_angle_scores[j].second;
//			++i_rslope;
//		}
//	}
//
//	vector<vector<int>> rslope_id_hierachy;
//	basicHierarchicalAlg<float>(rslope_angles, M_PI/180., rslope_id_hierachy);
//	vector<int>& first_class_id = rslope_id_hierachy[0];

	vector<float> angles_vector;
	vector<int> hist;
	vector<int> hist_smoothed;
	vector<vector<int>> hist_hierarchy;
	vector< pair<int, int> > id_regionid_lineid;
	int count_bins = 36;

	int temp_id = 0;
	for (int i = 0; i < count_char_regions; ++i) {
		CharRegion& region = charRegionArray[i];
		for (int j = 0; j < region.line_angles.size(); ++j) {
			angles_vector.push_back(region.line_angles[j]);
			id_regionid_lineid.push_back( make_pair(i, j) );
			temp_id++;
		}
	}

	pair<float, float> angle_range = make_pair(0.f, static_cast<float>(M_PI));
	clac1DHistogram<float>(angles_vector, count_bins, angle_range, hist,
			hist_hierarchy);
	conv1DSignal<int>(hist, {0.5f, 1.f, 0.5f}, hist_smoothed);

	for (int i = 0; i < count_bins; ++i) {
		cout << i << "--" << hist[i] << ": " << hist_smoothed[i];
		cout << endl;
	}

	// mean slope of this paper
	vector<int>::iterator it_max = std::max_element(hist_smoothed.begin(), hist_smoothed.end()) ;
	int slope_degree = std::distance(hist_smoothed.begin(), it_max);
	float mean_radians = M_PI/2.f + M_PI*float(slope_degree)/float(count_bins);
	float mean_slope = tan(mean_radians);
	if(mean_radians>0.f) {
		mean_slope = -mean_slope;
	}
	// get mean slope line for identification slope = rise/move = delta(y)/delta(x)
	int page_move = lines_im.cols;
	int page_rise = mean_slope * page_move;
	int x1_mean_line = 0;
	int x2_mean_line = page_move-1;
	int y1_mean_line = -page_rise/2 + lines_im.rows/2;
	int y2_mean_line = page_rise/2 + lines_im.rows/2;
	int x1_aux_line, x2_aux_line, y1_aux_line, y2_aux_line;

	line(lines_im, Point(x1_mean_line, y1_mean_line), Point(x2_mean_line, y2_mean_line),
							CharRegion::solarized_palette["blue"], 3);

	if(mean_radians<0.f) {
		// to avoid judgement of the direction of the point to the line we use this
		x1_aux_line = 0;
		y1_aux_line = 0;
		x2_aux_line = page_move;
		y2_aux_line = -page_rise;
	}
	else {
		x1_aux_line = 0;
		y1_aux_line = -page_rise;
		x2_aux_line = page_move;
		y2_aux_line = 0;
	}
	// pre-detemined base line
//	mrpt::math::TLine2D base_line_pre_det(
//											mrpt::math::TPoint2D(double(x1_aux_line), double(y1_aux_line)),
//											mrpt::math::TPoint2D(double(x2_aux_line), double(y2_aux_line))
//										);
	double A_aux_line, B_aux_line, C_aux_line;
	CharRegion::twoPoints_ABC(x1_aux_line, y1_aux_line, x2_aux_line, y2_aux_line,
					A_aux_line, B_aux_line, C_aux_line);

	//
	vector<double> distance_array;
	vector<int>& cand_line_ids = hist_hierarchy[slope_degree];
	for(auto& id: cand_line_ids) {
		pair<int, int>& regionid_lineid = id_regionid_lineid[id];
		int regionid = regionid_lineid.first;
		int lineid = regionid_lineid.second;
		CharRegion& cr = charRegionArray[regionid];

		int x_mean = (cr.bb_intersection_pts[lineid*4+0]+cr.bb_intersection_pts[lineid*4+2])/2;
		int y_mean = (cr.bb_intersection_pts[lineid*4+1]+cr.bb_intersection_pts[lineid*4+3])/2;

		line(lines_im,
				Point(cr.bb_intersection_pts[lineid*4+0], cr.bb_intersection_pts[lineid*4+1]),
				Point(cr.bb_intersection_pts[lineid*4+2], cr.bb_intersection_pts[lineid*4+3]),
									CharRegion::solarized_palette["blue"], 2);
		// mrpt::math::TPoint2D mean_pt(double(x_mean), double(y_mean));
		// double dist = base_line_pre_det.distance(mean_pt);	// fuck software development, i am not an expert
																// so I dont know why the mismatch from MRPT ducument
																// http://reference.mrpt.org/stable/structmrpt_1_1math_1_1_t_line2_d.html#a1484bc92c6516f9373049c5900d2f3a9
//		double a = base_line_pre_det.coefs[0];
//		double b = base_line_pre_det.coefs[1];
//		double c = base_line_pre_det.coefs[2];


		////  AFTER ALL, I USE MY OWN IMPLEMENTATION OF HIGH SCHOOL GEOMETRY
		double dist = abs( A_aux_line*x_mean+B_aux_line*y_mean+C_aux_line )/
							sqrt(A_aux_line*A_aux_line+B_aux_line*B_aux_line);
		// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

		distance_array.push_back(dist);	// EXTRA-COMPUTATION: vector
	}
	double sum_dist = 0.;
	std::for_each(distance_array.begin(), distance_array.end(),
			[&] (double n) {
				sum_dist += n;
			}
	);
	double mean_dist = sum_dist/distance_array.size();
	cout << "mean_dist : " << mean_dist << endl;

	float abs_slope_angle = abs(mean_slope);
	int descent = mean_dist / cos(abs_slope_angle);
	int x1_base_line = x1_aux_line;
	int x2_base_line = x2_aux_line;
	int y1_base_line = y1_aux_line+descent;
	int y2_base_line = y2_aux_line+descent;

	line(lines_im,
			Point(x1_base_line, y1_base_line),
			Point(x2_base_line, y2_base_line),
			CharRegion::solarized_palette["magenta"],
			5);

	double A_base_line, B_base_line, C_base_line;
	CharRegion::twoPoints_ABC(x1_base_line, y1_base_line,
								x2_base_line, y2_base_line,
								A_base_line, B_base_line, C_base_line);

	// get vertical line of base_line
	int mid_base_line_x = (x1_base_line + x2_base_line)/2;
	int mid_base_line_y = (y1_base_line + y2_base_line)/2;
	double A_base_line_vert, B_base_line_vert, C_base_line_vert;
	CharRegion::getVerticalLine(A_base_line, B_base_line, C_base_line,
									mid_base_line_x, mid_base_line_y,
									A_base_line_vert, B_base_line_vert, C_base_line_vert);

	cout << mid_base_line_x << " , " << mid_base_line_y << endl;
	int x1_base_line_vert, x2_base_line_vert, y1_base_line_vert, y2_base_line_vert;
	CharRegion::ABC_2points(A_base_line_vert, B_base_line_vert, C_base_line_vert,
					x1_base_line_vert, y1_base_line_vert, x2_base_line_vert, y2_base_line_vert);

	cout << A_base_line << "x + " << B_base_line  << "y + " << C_base_line
				<< " = 0" << endl;
	cout << x1_base_line << " , ";
	cout << y1_base_line << " , ";
	cout << x2_base_line << " , ";
	cout << y2_base_line << endl;
	cout << A_base_line_vert << "x + " << B_base_line_vert  << "y + " << C_base_line_vert
			<< " = 0" << endl;
	cout << x1_base_line_vert << " , ";
	cout << y1_base_line_vert << " , ";
	cout << x2_base_line_vert << " , ";
	cout << y2_base_line_vert << endl;

	line(lines_im,
			Point(x1_base_line_vert, y1_base_line_vert),
			Point(x2_base_line_vert, y2_base_line_vert),
			CharRegion::solarized_palette["orange"],
			5);

	imwrite("mean_lines.png", lines_im);

	// PageDeutschland page;
	PageDeutschland(A_base_line, B_base_line, C_base_line,
					A_base_line_vert, B_base_line_vert, C_base_line_vert,
					input.cols, input.rows,
					charRegionArray);

	//// free pres_mat
	delete[] pres_mat;

	// cout << "THE FUCK CRUSHED!?" << endl;

	return edge_swt;
}

void strokeWidthTransform(const Mat& edgeImage, Mat& gradientX, Mat& gradientY,
		bool dark_on_light, Mat& SWTImage, std::vector<Ray> & rays) {
	// First pass
	float prec = .05;
	for (int row = 0; row < edgeImage.rows; row++) {
		const uchar* ptr = (const uchar*) edgeImage.ptr(row);
		for (int col = 0; col < edgeImage.cols; col++) {
			if (*ptr > 0) {
				Ray r;

				SWTPoint2d p;
				p.x = col;
				p.y = row;
				r.p = p;
				std::vector<SWTPoint2d> points;
				points.push_back(p);

				float curX = (float) col + 0.5;
				float curY = (float) row + 0.5;
				int curPixX = col;
				int curPixY = row;
				float G_x = gradientX.at<float>(row, col);
				float G_y = gradientY.at<float>(row, col);
				// normalize gradient
				float mag = sqrt((G_x * G_x) + (G_y * G_y));
				if (dark_on_light) {
					G_x = -G_x / mag;
					G_y = -G_y / mag;
				} else {
					G_x = G_x / mag;
					G_y = G_y / mag;

				}
				while (true) {
					curX += G_x * prec;
					curY += G_y * prec;
					if ((int) (floor(curX)) != curPixX
							|| (int) (floor(curY)) != curPixY) {
						curPixX = (int) (floor(curX));
						curPixY = (int) (floor(curY));
						// check if pixel is outside boundary of image
						if (curPixX < 0 || (curPixX >= SWTImage.cols)
								|| curPixY < 0 || (curPixY >= SWTImage.rows)) {
							break;
						}
						SWTPoint2d pnew;
						pnew.x = curPixX;
						pnew.y = curPixY;
						points.push_back(pnew);

						if (edgeImage.at < uchar > (curPixY, curPixX) > 0) {
							r.q = pnew;
							// dot product
							float G_xt = gradientX.at<float>(curPixY, curPixX);
							float G_yt = gradientY.at<float>(curPixY, curPixX);
							mag = sqrt((G_xt * G_xt) + (G_yt * G_yt));
							if (dark_on_light) {
								G_xt = -G_xt / mag;
								G_yt = -G_yt / mag;
							} else {
								G_xt = G_xt / mag;
								G_yt = G_yt / mag;

							}

							if (acos(G_x * -G_xt + G_y * -G_yt) < PI / 2.0) {
								float length =
										sqrt(
												((float) r.q.x - (float) r.p.x)
														* ((float) r.q.x
																- (float) r.p.x)
														+ ((float) r.q.y
																- (float) r.p.y)
																* ((float) r.q.y
																		- (float) r.p.y));
								for (std::vector<SWTPoint2d>::iterator pit =
										points.begin(); pit != points.end();
										pit++) {
									if (SWTImage.at<float>(pit->y, pit->x)
											< 0) {
										SWTImage.at<float>(pit->y, pit->x) =
												length;
									} else {
										SWTImage.at<float>(pit->y, pit->x) =
												std::min(length,
														SWTImage.at<float>(
																pit->y,
																pit->x));
									}
								}
								r.points = points;
								rays.push_back(r);
							}
							break;
						}
					}
				}
			}
			ptr++;
		}
	}

}

void SWTMedianFilter(Mat& SWTImage, std::vector<Ray> & rays) {
	for (auto& rit : rays) {
		for (auto& pit : rit.points) {
			pit.SWT = SWTImage.at<float>(pit.y, pit.x);
		}
		std::sort(rit.points.begin(), rit.points.end(), &Point2dSort);
		float median = (rit.points[rit.points.size() / 2]).SWT;
		for (auto& pit : rit.points) {
			SWTImage.at<float>(pit.y, pit.x) = std::min(pit.SWT, median);
		}
	}
}

bool Point2dSort(const SWTPoint2d &lhs, const SWTPoint2d &rhs) {
	return lhs.SWT < rhs.SWT;
}

std::vector<std::vector<SWTPoint2d> > findLegallyConnectedComponents(
		Mat& SWTImage, std::vector<Ray> & rays) {
	boost::unordered_map<int, int> map;
	boost::unordered_map<int, SWTPoint2d> revmap;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	int num_vertices = 0;
	// Number vertices for graph.  Associate each point with number
	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*) SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				map[row * SWTImage.cols + col] = num_vertices;
				SWTPoint2d p;
				p.x = col;
				p.y = row;
				revmap[num_vertices] = p;
				num_vertices++;
			}
			ptr++;
		}
	}

	Graph g(num_vertices);

	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*) SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				// check pixel to the right, right-down, down, left-down
				int this_pixel = map[row * SWTImage.cols + col];
				if (col + 1 < SWTImage.cols) {
					float right = SWTImage.at<float>(row, col + 1);
					if (right > 0
							&& ((*ptr) / right <= 3.0 || right / (*ptr) <= 3.0))
						boost::add_edge(this_pixel,
								map.at(row * SWTImage.cols + col + 1), g);
				}
				if (row + 1 < SWTImage.rows) {
					if (col + 1 < SWTImage.cols) {
						float right_down = SWTImage.at<float>(row + 1, col + 1);
						if (right_down > 0
								&& ((*ptr) / right_down <= 3.0
										|| right_down / (*ptr) <= 3.0))
							boost::add_edge(this_pixel,
									map.at((row + 1) * SWTImage.cols + col + 1),
									g);
					}
					float down = SWTImage.at<float>(row + 1, col);
					if (down > 0
							&& ((*ptr) / down <= 3.0 || down / (*ptr) <= 3.0))
						boost::add_edge(this_pixel,
								map.at((row + 1) * SWTImage.cols + col), g);
					if (col - 1 >= 0) {
						float left_down = SWTImage.at<float>(row + 1, col - 1);
						if (left_down > 0
								&& ((*ptr) / left_down <= 3.0
										|| left_down / (*ptr) <= 3.0))
							boost::add_edge(this_pixel,
									map.at((row + 1) * SWTImage.cols + col - 1),
									g);
					}
				}
			}
			ptr++;
		}
	}

	std::vector<int> c(num_vertices);

	int num_comp = connected_components(g, &c[0]);

	std::vector < std::vector<SWTPoint2d> > components;
	components.reserve(num_comp);
	std::cout << "Before filtering, " << num_comp << " components and "
			<< num_vertices << " vertices" << std::endl;
	for (int j = 0; j < num_comp; j++) {
		std::vector<SWTPoint2d> tmp;
		components.push_back(tmp);
	}
	for (int j = 0; j < num_vertices; j++) {
		SWTPoint2d p = revmap[j];
		(components[c[j]]).push_back(p);
	}

	return components;
}

std::vector<std::vector<SWTPoint2d> > findLegallyConnectedComponentsRAY(
		Mat& SWTImage, std::vector<Ray> & rays) {
	boost::unordered_map<int, int> map;
	boost::unordered_map<int, SWTPoint2d> revmap;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	int num_vertices = 0;
	// Number vertices for graph.  Associate each point with number
	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*) SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				map[row * SWTImage.cols + col] = num_vertices;
				SWTPoint2d p;
				p.x = col;
				p.y = row;
				revmap[num_vertices] = p;
				num_vertices++;
			}
			ptr++;
		}
	}

	Graph g(num_vertices);

	// Traverse and add edges to graph
	for (std::vector<Ray>::const_iterator it = rays.begin(); it != rays.end();
			it++) {
		float lastSW = 0;
		int lastRow = 0;
		int lastCol = 0;
		for (std::vector<SWTPoint2d>::const_iterator it2 = it->points.begin();
				it2 != it->points.end(); it2++) {
			float currentSW = SWTImage.at<float>(it2->y, it2->x);
			if (lastSW == 0) {
			} else if (lastSW / currentSW <= 3.0 || currentSW / lastSW <= 3.0) {
				boost::add_edge(map.at(it2->y * SWTImage.cols + it2->x),
						map.at(lastRow * SWTImage.cols + lastCol), g);
			}
			lastSW = currentSW;
			lastRow = it2->y;
			lastCol = it2->x;
		}
		lastSW = 0;
		lastRow = 0;
		lastCol = 0;
	}

	std::vector<int> c(num_vertices);

	int num_comp = connected_components(g, &c[0]);

	std::vector < std::vector<SWTPoint2d> > components;
	components.reserve(num_comp);
	std::cout << "Before filtering, " << num_comp << " components and "
			<< num_vertices << " vertices" << std::endl;
	for (int j = 0; j < num_comp; j++) {
		std::vector<SWTPoint2d> tmp;
		components.push_back(tmp);
	}
	for (int j = 0; j < num_vertices; j++) {
		SWTPoint2d p = revmap[j];
		(components[c[j]]).push_back(p);
	}

	return components;
}

void componentStats(Mat& SWTImage, const std::vector<SWTPoint2d> & component,
		float & mean, float & variance, float & median, int & minx, int & miny,
		int & maxx, int & maxy) {
	std::vector<float> temp;
	temp.reserve(component.size());
	mean = 0;
	variance = 0;
	minx = 1000000;
	miny = 1000000;
	maxx = 0;
	maxy = 0;
	for (std::vector<SWTPoint2d>::const_iterator it = component.begin();
			it != component.end(); it++) {
		float t = SWTImage.at<float>(it->y, it->x);
		mean += t;
		temp.push_back(t);
		miny = std::min(miny, it->y);
		minx = std::min(minx, it->x);
		maxy = std::max(maxy, it->y);
		maxx = std::max(maxx, it->x);
	}
	mean = mean / ((float) component.size());
	for (std::vector<float>::const_iterator it = temp.begin(); it != temp.end();
			it++) {
		variance += (*it - mean) * (*it - mean);
	}
	variance = variance / ((float) component.size());
	std::sort(temp.begin(), temp.end());
	median = temp[temp.size() / 2];
}

void filterComponents(Mat& SWTImage,
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<std::vector<SWTPoint2d> > & validComponents,
		std::vector<Point2dFloat> & compCenters,
		std::vector<float> & compMedians,
		std::vector<SWTPoint2d> & compDimensions,
		std::vector<SWTPointPair2d> & compBB) {
	validComponents.reserve(components.size());
	compCenters.reserve(components.size());
	compMedians.reserve(components.size());
	compDimensions.reserve(components.size());
	// bounding boxes
	compBB.reserve(components.size());
	for (std::vector<std::vector<SWTPoint2d> >::iterator it =
			components.begin(); it != components.end(); it++) {
		// compute the stroke width mean, variance, median
		float mean, variance, median;
		int minx, miny, maxx, maxy;
		componentStats(SWTImage, (*it), mean, variance, median, minx, miny,
				maxx, maxy);

		// check if variance is less than half the mean
		if (variance > 0.5 * mean) {
			continue;
		}

		float length = (float) (maxx - minx + 1);
		float width = (float) (maxy - miny + 1);

		// check font height
		if (width > 300) {
			continue;
		}

		float area = length * width;
		float rminx = (float) minx;
		float rmaxx = (float) maxx;
		float rminy = (float) miny;
		float rmaxy = (float) maxy;
		// compute the rotated bounding box
		float increment = 1. / 36.;
		for (float theta = increment * PI; theta < PI / 2.0;
				theta += increment * PI) {
			float xmin, xmax, ymin, ymax, xtemp, ytemp, ltemp, wtemp;
			xmin = 1000000;
			ymin = 1000000;
			xmax = 0;
			ymax = 0;
			for (unsigned int i = 0; i < (*it).size(); i++) {
				xtemp = (*it)[i].x * cos(theta) + (*it)[i].y * -sin(theta);
				ytemp = (*it)[i].x * sin(theta) + (*it)[i].y * cos(theta);
				xmin = std::min(xtemp, xmin);
				xmax = std::max(xtemp, xmax);
				ymin = std::min(ytemp, ymin);
				ymax = std::max(ytemp, ymax);
			}
			ltemp = xmax - xmin + 1;
			wtemp = ymax - ymin + 1;
			if (ltemp * wtemp < area) {
				area = ltemp * wtemp;
				length = ltemp;
				width = wtemp;
			}
		}
		// check if the aspect ratio is between 1/10 and 10
		if (length / width < 1. / 10. || length / width > 10.) {
			continue;
		}

		// compute the diameter TODO finish
		// compute dense representation of component
		std::vector < std::vector<float> > denseRepr;
		denseRepr.reserve(maxx - minx + 1);
		for (int i = 0; i < maxx - minx + 1; i++) {
			std::vector<float> tmp;
			tmp.reserve(maxy - miny + 1);
			denseRepr.push_back(tmp);
			for (int j = 0; j < maxy - miny + 1; j++) {
				denseRepr[i].push_back(0);
			}
		}
		for (std::vector<SWTPoint2d>::iterator pit = it->begin();
				pit != it->end(); pit++) {
			(denseRepr[pit->x - minx])[pit->y - miny] = 1;
		}
		// create graph representing components
		const int num_nodes = it->size();
		/*
		 E edges[] = { E(0,2),
		 E(1,1), E(1,3), E(1,4),
		 E(2,1), E(2,3),
		 E(3,4),
		 E(4,0), E(4,1) };

		 Graph G(edges + sizeof(edges) / sizeof(E), weights, num_nodes);
		 */
		Point2dFloat center;
		center.x = ((float) (maxx + minx)) / 2.0;
		center.y = ((float) (maxy + miny)) / 2.0;

		SWTPoint2d dimensions;
		dimensions.x = maxx - minx + 1;
		dimensions.y = maxy - miny + 1;

		SWTPoint2d bb1;
		bb1.x = minx;
		bb1.y = miny;

		SWTPoint2d bb2;
		bb2.x = maxx;
		bb2.y = maxy;
		SWTPointPair2d pair(bb1, bb2);

		compBB.push_back(pair);
		compDimensions.push_back(dimensions);
		compMedians.push_back(median);
		compCenters.push_back(center);
		validComponents.push_back(*it);
	}
	std::vector < std::vector<SWTPoint2d> > tempComp;
	std::vector<SWTPoint2d> tempDim;
	std::vector<float> tempMed;
	std::vector<Point2dFloat> tempCenters;
	std::vector<SWTPointPair2d> tempBB;
	tempComp.reserve(validComponents.size());
	tempCenters.reserve(validComponents.size());
	tempDim.reserve(validComponents.size());
	tempMed.reserve(validComponents.size());
	tempBB.reserve(validComponents.size());
	for (unsigned int i = 0; i < validComponents.size(); i++) {
		int count = 0;
		for (unsigned int j = 0; j < validComponents.size(); j++) {
			if (i != j) {
				if (compBB[i].first.x <= compCenters[j].x
						&& compBB[i].second.x >= compCenters[j].x
						&& compBB[i].first.y <= compCenters[j].y
						&& compBB[i].second.y >= compCenters[j].y) {
					count++;
				}
			}
		}
		if (count < 2) {
			tempComp.push_back(validComponents[i]);
			tempCenters.push_back(compCenters[i]);
			tempMed.push_back(compMedians[i]);
			tempDim.push_back(compDimensions[i]);
			tempBB.push_back(compBB[i]);
		}
	}
	validComponents = tempComp;
	compDimensions = tempDim;
	compMedians = tempMed;
	compCenters = tempCenters;
	compBB = tempBB;

	compDimensions.reserve(tempComp.size());
	compMedians.reserve(tempComp.size());
	compCenters.reserve(tempComp.size());
	validComponents.reserve(tempComp.size());
	compBB.reserve(tempComp.size());

	std::cout << "After filtering " << validComponents.size() << " components"
			<< std::endl;
}

bool sharesOneEnd(Chain c0, Chain c1) {
	if (c0.p == c1.p || c0.p == c1.q || c0.q == c1.q || c0.q == c1.p) {
		return true;
	} else {
		return false;
	}
}

bool chainSortDist(const Chain &lhs, const Chain &rhs) {
	return lhs.dist < rhs.dist;
}

bool chainSortLength(const Chain &lhs, const Chain &rhs) {
	return lhs.components.size() > rhs.components.size();
}

std::vector<Chain> makeChains(const Mat& colorImage,
		std::vector<std::vector<SWTPoint2d> > & components,
		std::vector<Point2dFloat> & compCenters,
		std::vector<float> & compMedians,
		std::vector<SWTPoint2d> & compDimensions,
		std::vector<SWTPointPair2d> & compBB) {
	assert(compCenters.size() == components.size());
	// make vector of color averages
	std::vector<Point3dFloat> colorAverages;
	colorAverages.reserve(components.size());
	for (std::vector<std::vector<SWTPoint2d> >::iterator it =
			components.begin(); it != components.end(); it++) {
		Point3dFloat mean;
		mean.x = 0;
		mean.y = 0;
		mean.z = 0;
		int num_points = 0;
		for (std::vector<SWTPoint2d>::iterator pit = it->begin();
				pit != it->end(); pit++) {
			mean.x += (float) colorImage.at < uchar > (pit->y, (pit->x) * 3);
			mean.y += (float) colorImage.at < uchar
					> (pit->y, (pit->x) * 3 + 1);
			mean.z += (float) colorImage.at < uchar
					> (pit->y, (pit->x) * 3 + 2);
			num_points++;
		}
		mean.x = mean.x / ((float) num_points);
		mean.y = mean.y / ((float) num_points);
		mean.z = mean.z / ((float) num_points);
		colorAverages.push_back(mean);
	}

	// form all eligible pairs and calculate the direction of each
	std::vector<Chain> chains;
	for (unsigned int i = 0; i < components.size(); i++) {
		for (unsigned int j = i + 1; j < components.size(); j++) {
			// TODO add color metric
			if ((compMedians[i] / compMedians[j] <= 2.0
					|| compMedians[j] / compMedians[i] <= 2.0)
					&& (compDimensions[i].y / compDimensions[j].y <= 2.0
							|| compDimensions[j].y / compDimensions[i].y <= 2.0)) {
				float dist = (compCenters[i].x - compCenters[j].x)
						* (compCenters[i].x - compCenters[j].x)
						+ (compCenters[i].y - compCenters[j].y)
								* (compCenters[i].y - compCenters[j].y);
				float colorDist = (colorAverages[i].x - colorAverages[j].x)
						* (colorAverages[i].x - colorAverages[j].x)
						+ (colorAverages[i].y - colorAverages[j].y)
								* (colorAverages[i].y - colorAverages[j].y)
						+ (colorAverages[i].z - colorAverages[j].z)
								* (colorAverages[i].z - colorAverages[j].z);
				if (dist
						< 9
								* (float) (std::max(
										std::min(compDimensions[i].x,
												compDimensions[i].y),
										std::min(compDimensions[j].x,
												compDimensions[j].y)))
								* (float) (std::max(
										std::min(compDimensions[i].x,
												compDimensions[i].y),
										std::min(compDimensions[j].x,
												compDimensions[j].y)))
						&& colorDist < 1600) {
					Chain c;
					c.p = i;
					c.q = j;
					std::vector<int> comps;
					comps.push_back(c.p);
					comps.push_back(c.q);
					c.components = comps;
					c.dist = dist;
					float d_x = (compCenters[i].x - compCenters[j].x);
					float d_y = (compCenters[i].y - compCenters[j].y);
					/*
					 float d_x = (compBB[i].first.x - compBB[j].second.x);
					 float d_y = (compBB[i].second.y - compBB[j].second.y);
					 */
					float mag = sqrt(d_x * d_x + d_y * d_y);
					d_x = d_x / mag;
					d_y = d_y / mag;
					Point2dFloat dir;
					dir.x = d_x;
					dir.y = d_y;
					c.direction = dir;
					chains.push_back(c);

					/*std::cerr << c.p << " " << c.q << std::endl;
					 std::cerr << c.direction.x << " " << c.direction.y << std::endl;
					 std::cerr << compCenters[c.p].x << " " << compCenters[c.p].y << std::endl;
					 std::cerr << compCenters[c.q].x << " " << compCenters[c.q].y << std::endl;
					 std::cerr << std::endl;
					 std::cerr << colorDist << std::endl; */
				}
			}
		}
	}
	std::cout << chains.size() << " eligible pairs" << std::endl;
	std::sort(chains.begin(), chains.end(), &chainSortDist);

	std::cerr << std::endl;
	const float strictness = PI / 6.0;
	//merge chains
	int merges = 1;
	while (merges > 0) {
		for (unsigned int i = 0; i < chains.size(); i++) {
			chains[i].merged = false;
		}
		merges = 0;
		std::vector<Chain> newchains;
		for (unsigned int i = 0; i < chains.size(); i++) {
			for (unsigned int j = 0; j < chains.size(); j++) {
				if (i != j) {
					if (!chains[i].merged && !chains[j].merged
							&& sharesOneEnd(chains[i], chains[j])) {
						if (chains[i].p == chains[j].p) {
							if (acos(
									chains[i].direction.x
											* -chains[j].direction.x
											+ chains[i].direction.y
													* -chains[j].direction.y)
									< strictness) {
								/*      if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 }
								 std::cerr << 1 <<std::endl;

								 std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								 std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								 std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								 std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								 std::cerr << compCenters[chains[j].q].x << " " << compCenters[chains[j].q].y << std::endl;
								 std::cerr << std::endl; */

								chains[i].p = chains[j].q;
								for (std::vector<int>::iterator it =
										chains[j].components.begin();
										it != chains[j].components.end();
										it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x
										- compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y
										- compCenters[chains[i].q].y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								float mag = sqrt(d_x * d_x + d_y * d_y);
								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								 i=0;
								 if (i == chains.size() - 1) i=-1;
								 std::stable_sort(chains.begin(), chains.end(), &chainSortLength);*/
							}
						} else if (chains[i].p == chains[j].q) {
							if (acos(
									chains[i].direction.x
											* chains[j].direction.x
											+ chains[i].direction.y
													* chains[j].direction.y)
									< strictness) {
								/*
								 if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 }
								 std::cerr << 2 <<std::endl;

								 std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								 std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								 std::cerr << chains[i].direction.x << " " << chains[i].direction.y << std::endl;
								 std::cerr << chains[j].direction.x << " " << chains[j].direction.y << std::endl;
								 std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								 std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								 std::cerr << compCenters[chains[j].p].x << " " << compCenters[chains[j].p].y << std::endl;
								 std::cerr << std::endl; */

								chains[i].p = chains[j].p;
								for (std::vector<int>::iterator it =
										chains[j].components.begin();
										it != chains[j].components.end();
										it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x
										- compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y
										- compCenters[chains[i].q].y);
								float mag = sqrt(d_x * d_x + d_y * d_y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								d_x = d_x / mag;
								d_y = d_y / mag;

								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								 i=0;
								 if (i == chains.size() - 1) i=-1;
								 std::stable_sort(chains.begin(), chains.end(), &chainSortLength); */
							}
						} else if (chains[i].q == chains[j].p) {
							if (acos(
									chains[i].direction.x
											* chains[j].direction.x
											+ chains[i].direction.y
													* chains[j].direction.y)
									< strictness) {
								/*                           if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 }
								 std::cerr << 3 <<std::endl;

								 std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								 std::cerr << chains[j].p << " " << chains[j].q << std::endl;

								 std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								 std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								 std::cerr << compCenters[chains[j].q].x << " " << compCenters[chains[j].q].y << std::endl;
								 std::cerr << std::endl; */
								chains[i].q = chains[j].q;
								for (std::vector<int>::iterator it =
										chains[j].components.begin();
										it != chains[j].components.end();
										it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x
										- compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y
										- compCenters[chains[i].q].y);
								float mag = sqrt(d_x * d_x + d_y * d_y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;

								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								 i=0;
								 if (i == chains.size() - 1) i=-1;
								 std::stable_sort(chains.begin(), chains.end(), &chainSortLength); */
							}
						} else if (chains[i].q == chains[j].q) {
							if (acos(
									chains[i].direction.x
											* -chains[j].direction.x
											+ chains[i].direction.y
													* -chains[j].direction.y)
									< strictness) {
								/*           if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 } else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								 std::cout << "CRAZY ERROR" << std::endl;
								 }
								 std::cerr << 4 <<std::endl;
								 std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								 std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								 std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								 std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								 std::cerr << compCenters[chains[j].p].x << " " << compCenters[chains[j].p].y << std::endl;
								 std::cerr << std::endl; */
								chains[i].q = chains[j].p;
								for (std::vector<int>::iterator it =
										chains[j].components.begin();
										it != chains[j].components.end();
										it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x
										- compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y
										- compCenters[chains[i].q].y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								float mag = sqrt(d_x * d_x + d_y * d_y);
								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								 i=0;
								 if (i == chains.size() - 1) i=-1;
								 std::stable_sort(chains.begin(), chains.end(), &chainSortLength);*/
							}
						}
					}
				}
			}
		}
		for (unsigned int i = 0; i < chains.size(); i++) {
			if (!chains[i].merged) {
				newchains.push_back(chains[i]);
			}
		}
		chains = newchains;
		std::stable_sort(chains.begin(), chains.end(), &chainSortLength);
	}

	std::vector<Chain> newchains;
	newchains.reserve(chains.size());
	for (std::vector<Chain>::iterator cit = chains.begin(); cit != chains.end();
			cit++) {
		if (cit->components.size() >= 3) {
			newchains.push_back(*cit);
		}
	}
	chains = newchains;
	std::cout << chains.size() << " chains after merging" << std::endl;
	return chains;
}

}
