#include <iostream>
#include <vector>

#include "CharRegion.h"
#include <opencv2/opencv.hpp>

#include <assert.h>

#ifndef PAGEPARTITION_H
#define PAGEPARTITION_H

using namespace std;
using namespace cv;

class PageDistrict {
public:
//	PageDistrict(vector<CharRegion>& charRegionArray) {
//		p_charRegionArray = &charRegionArray;
//	}
  vector<pair<int, int>> lines_array;

  vector<pair<int, vector<double>>> id_data_array;      // [id, {dist_value, angle, x, y}]
  vector<pair<int, vector<double>>> id_data_array2;// [id, {dist_value, angle, x, y}]
  vector<pair<int, vector<double>>> id_data_array3;// [id, {dist_value, angle, x, y}]
  vector<CharRegion>* p_charRegionArray;
  vector<bool> valid_table;

  double A_base_line, B_base_line, C_base_line;
  double A_base_line_vert, B_base_line_vert, C_base_line_vert;
  int x_center, y_center;

public:
  PageDistrict();

  void setRegionArray(vector<CharRegion>& charRegionArray, double g_A_base_line,
      double g_B_base_line, double g_C_base_line, double g_A_base_line_vert,
      double g_B_base_line_vert, double g_C_base_line_vert);

  void getSlopeArray(double* north_unit_vector);
  void bullshitsFilter();
  void oddDistanceFilter(int scanRange,
      vector<pair<int, vector<double>>>& source_signal,
      double odd_angle_thresh=6.,
      int count_odd_thresh=3);

  bool tellIfSheered(double mean_sheered_angle=3.*M_PI/180);
};

class PageDeutschland {
public:

  PageDeutschland(double g_A_base_line, double g_B_base_line, double g_C_base_line,
      double g_A_base_line_vert, double g_B_base_line_vert, double g_C_base_line_vert,
      int g_ori_width, int g_ori_height, vector<CharRegion>& charRegionArray);

  void getWolfsburgRegularPts(const Mat& im_in);
  void getJenaRegularPts(const Mat& im_in);
  void getStuttgartRegularPts(const Mat& im_in);
  void getMunchenRegularPts(const Mat& im_in);


public:
  PageDistrict regions[4];
  PageDistrict *p_wolfsburg, *p_jena, *p_stuttgart, *p_munchen;
  int code_wolfsburg, code_jena, code_stuttgart, code_munchen;
  double A_base_line, B_base_line, C_base_line;
  double A_base_line_vert, B_base_line_vert, C_base_line_vert;
  int ori_width, ori_height;

  double x1_base_line, y1_base_line, x2_base_line, y2_base_line;
  double x1_base_line_vert, y1_base_line_vert, x2_base_line_vert, y2_base_line_vert;
  int x_base_midpoint, y_base_midpoint;

  double A_base_line_west, B_base_line_west, C_base_line_west;
  double A_base_line_east, B_base_line_east, C_base_line_east;
  double A_base_line_vert_north, B_base_line_vert_north, C_base_line_vert_north;
  double A_base_line_vert_south, B_base_line_vert_south, C_base_line_vert_south;

  double north_unit_vector[2];

  template<typename Dtype>
  bool isWest(Dtype x, Dtype y);

  template<typename Dtype>
  bool isEast(Dtype x, Dtype y);

  template<typename Dtype>
  bool isNorth(Dtype x, Dtype y);

  template<typename Dtype>
  bool isSouth(Dtype x, Dtype y);

  char judgeNPE(const double& A, const double& B, const double& C, const double& x,
      const double& y);

  pair<double, double> stern_des_sudens;
  pair<double, double> stern_des_nordens;
  pair<double, double> stern_des_ostens;
  pair<double, double> stern_des_westens;

  pair<double, double> star_of_northwest;     // maybe exceed the range of image pixel
  pair<double, double> star_of_northeast;     // maybe exceed the range of image pixel
  pair<double, double> star_of_southwest;     // maybe exceed the range of image pixel
  pair<double, double> star_of_southeast;     // maybe exceed the range of image pixel

  Mat wolfsburg_im;
  Mat jena_im;
  Mat stuttgart_im;
  Mat munchen_im;

private:
  void genDeDistrictCode(int ** corner_pts);

  void reorganizeRegions(vector<CharRegion>& charRegionArray);

  // return: low 4 bits: 0x0 --> NPE1==N; 0x1 --> NPE1==P
  //			high 4 bits: 0x0 --> NPE2==N; 0x1 --> NPE2==P
  int encode2NPE(char NPE1, char NPE2);

  char north_code;
  char south_code;
  char west_code;
  char east_code;

};

#endif
