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
#include <math.h>
#include <stdlib.h>

#include <unordered_map>
#include <utility>

// #define _USE_MATH_DEFINES

using namespace std;

class CharRegion {
public:
	CharRegion(const int arg_count_pixs, const int arg_count_contours,
			const int arg_count_edges, const int arg_x_max, const int arg_x_min,
			const int arg_y_max, const int arg_y_min);

	void push_region_pix(const int arg_x, const int arg_y);

	void push_region_contours(const int arg_x, const int arg_y);

	void push_edges(const int arg_x, const int arg_y);

	void plot_char_region(const string& filename);

	bool is_measurable();

	void ransac_find_lines();

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

	int const_nlines;
	vector<double> line_paramsABC; // A, B, C
	vector<int> line_paramsPts;	// X1, Y1, X2, Y2
	vector<double> line_paramsABC_global;	// A_global, B_global, C_global
	vector<int> line_paramsPts_global;		// X1_global, Y1_global,
	vector<int> bb_intersection_pts;	// x1, y1, x2, y2, ...
	vector<int> line_pts;

	vector<float> line_slopes;
	vector<float> line_angles;

	vector<vector<int>> line_clusters_ids;

//	vector<float> region_angles;		// the following 2 are index draw from line_clusters_ids
//	vector<float> region_angles_score;

	vector<pair<float, float>> region_angle_scores;

	static unordered_map<string, cv::Scalar> solarized_palette;

	template <typename T>
	static T get_angle_diff(T a1, T a2);

private:
	// auxiliary var
	int icap_pixs;
	int icap_contours;
	int icap_edges;

	static void ABC_2points(const double A, const double B, const double C,
			int& x1, int& y1, int& x2, int& y2);

	static void ABC_2points(const double A, const double B, const double C,
			double& x1, double& y1, double& x2, double& y2);

	void cvtLocal2Global();

	void calcIntersectionLine();

	void refine_lines();

	void get_slope_angle(double A, double B, double C, float& slope,
			float& angle);

	// float get_angle_diff(float a1, float a2);

	void calc_lines_info();
};

template<typename T>
void partition1dData(const vector<T>& Data, const T distance_thresh,
		vector<vector<int>>& id_hierachy);

template<typename T>
void basicalHierarchicalAlg(const vector<T>& data, const T diff_thresh,
		vector<vector<int>>& id_hierachy);

#endif /* INCLUDE_CHARREGION_H_ */
