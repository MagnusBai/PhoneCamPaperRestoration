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
  vector<pair<int, vector<double>>> id_data_array2;      // [id, {dist_value, angle, x, y}]
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
};

class PageDeutschland {
public:

  PageDeutschland(double g_A_base_line, double g_B_base_line, double g_C_base_line,
      double g_A_base_line_vert, double g_B_base_line_vert, double g_C_base_line_vert,
      int g_ori_width, int g_ori_height, vector<CharRegion>& charRegionArray);

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

  double A_base_line_left, B_base_line_left, C_base_line_left;
  double A_base_line_right, B_base_line_right, C_base_line_right;

  double north_unit_vector[2];

  char judgeNPE(const double& A, const double& B, const double& C, const double& x,
      const double& y);

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

};

#endif