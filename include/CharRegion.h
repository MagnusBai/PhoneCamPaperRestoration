/*
 * CharRegion.h
 *
 *  Created on: Dec 18, 2016
 *      Author: magnus
 */

#ifndef INCLUDE_CHARREGION_H_
#define INCLUDE_CHARREGION_H_

#include <vector>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <string>

using namespace std;

class CharRegion {
public:
	CharRegion(const int arg_count_pixs, const int arg_count_contours,
			const int arg_count_edges, const int arg_x_max, const int arg_x_min,
			const int arg_y_max, const int arg_y_min) :
			count_pixs(arg_count_pixs), count_contours(arg_count_contours), count_edges(
					arg_count_edges), x_max(arg_x_max), x_min(arg_x_min), y_max(
					arg_y_max), y_min(arg_y_min), icap_pixs(0), icap_contours(
					0), icap_edges(0) {

		region_pixs.resize(count_pixs * 2, 0);
		region_contours.resize(count_contours * 2, 0);
		char_edge.resize(count_edges * 2, 0);

	}

	void push_region_pix(const int arg_x, const int arg_y) {
		region_pixs[2 * icap_pixs] = arg_x;
		region_pixs[2 * icap_pixs + 1] = arg_y;
		++icap_pixs;
	}

	void push_region_contours(const int arg_x, const int arg_y) {
		region_contours[2 * icap_contours] = arg_x;
		region_contours[2 * icap_contours + 1] = arg_y;
		++icap_contours;
	}

	void push_edges(const int arg_x, const int arg_y) {
		char_edge[2 * icap_edges] = arg_x;
		char_edge[2 * icap_edges + 1] = arg_y;
		++icap_edges;
	}

	void plot_char_region(const string& filename) {
		cout << y_max << " , " << y_min << endl;
		cout << x_max << " , " << x_min << endl;
		cv::Mat im(y_max - y_min + 1, x_max - x_min + 1, CV_8UC1,
				cv::Scalar(0));

		for (int i = 0; i < count_pixs; ++i) {
			int x = region_pixs[i * 2] - x_min;
			int y = region_pixs[i * 2 + 1] - y_min;
			im.at < uchar > (y, x) = 150;
		}
//		for (int i = 0; i < count_contours; ++i) {
//			int x = region_contours[i * 2] - x_min;
//			int y = region_contours[i * 2 + 1] - y_min;
//			im.at < uchar > (y, x) = 255;
//		}
		for (int i = 0; i < count_edges; ++i) {
			int x = char_edge[i * 2] - x_min;
			int y = char_edge[i * 2 + 1] - y_min;
			im.at < uchar > (y, x) = 255;
		}

		cv::imwrite(filename.c_str(), im);
	}

	void log();

	// data
	vector<int> region_pixs;
	vector<int> region_contours;
	vector<int> char_edge;
	const int count_pixs;
	const int count_contours;
	const int count_edges;
	const int x_max;
	const int x_min;
	const int y_max;
	const int y_min;

	// auxiliary var
	int icap_pixs;
	int icap_contours;
	int icap_edges;
};

#endif /* INCLUDE_CHARREGION_H_ */
