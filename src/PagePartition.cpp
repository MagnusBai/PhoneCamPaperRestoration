/*
 * PagePartition.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: magnus
 */

#include <algorithm>

#include "PagePartition.h"

PageDeutschland::PageDeutschland(double g_A_base_line, double g_B_base_line, double g_C_base_line,
    double g_A_base_line_vert, double g_B_base_line_vert, double g_C_base_line_vert, int g_ori_width, int g_ori_height,
    vector<CharRegion>& charRegionArray) :
    A_base_line(g_A_base_line), B_base_line(g_B_base_line), C_base_line(g_C_base_line), A_base_line_vert(
        g_A_base_line_vert), B_base_line_vert(g_B_base_line_vert), C_base_line_vert(g_C_base_line_vert), p_wolfsburg(
        &regions[0]), p_jena(&regions[1]), p_stuttgart(&regions[2]), p_munchen(&regions[3]), ori_width(g_ori_width), ori_height(
        g_ori_height) {
  // get info of base line
  CharRegion::ABC_2points(A_base_line, B_base_line, C_base_line, x1_base_line, y1_base_line, x2_base_line, y2_base_line,
      0, g_ori_width - 1, 0, g_ori_height - 1);
  CharRegion::ABC_2points(A_base_line_vert, B_base_line_vert, C_base_line_vert, x1_base_line_vert, y1_base_line_vert,
      x2_base_line_vert, y2_base_line_vert, 0, g_ori_width - 1, 0, g_ori_height - 1);
  bool temp_flag;
  CharRegion::getIntersectionPoint(A_base_line, B_base_line, C_base_line, A_base_line_vert, B_base_line_vert,
      C_base_line_vert, x_base_midpoint, y_base_midpoint, temp_flag);

  double north_delta_x = x_base_midpoint - x1_base_line_vert;
  double north_delta_y = y_base_midpoint - y1_base_line_vert;
  double north_base = sqrt(north_delta_x * north_delta_x + north_delta_y * north_delta_y);
  north_unit_vector[0] = north_delta_x / north_base;
  north_unit_vector[1] = north_delta_y / north_base;

  p_wolfsburg->setRegionArray(charRegionArray, A_base_line, B_base_line, C_base_line, A_base_line_vert,
      B_base_line_vert, C_base_line_vert);
  p_jena->setRegionArray(charRegionArray, A_base_line, B_base_line, C_base_line, A_base_line_vert, B_base_line_vert,
      C_base_line_vert);
  p_stuttgart->setRegionArray(charRegionArray, A_base_line, B_base_line, C_base_line, A_base_line_vert,
      B_base_line_vert, C_base_line_vert);
  p_munchen->setRegionArray(charRegionArray, A_base_line, B_base_line, C_base_line, A_base_line_vert, B_base_line_vert,
      C_base_line_vert);

  // 4 corner
  // Wolfsburg (north-west)
  int wolfsburg_pts[6] = { 0, 0, -1, 0, 0, -1 };
  // Jena (north-east)
  int jena_pts[6] = { ori_width, 0, ori_width + 1, 0, ori_width, -1 };
  // Stuttgart (south-west)
  int stuttgart_pts[6] = { 0, ori_height, 0, ori_height + 1, -1, ori_height };
  // München (south-east)
  int muchen_pts[6] = { ori_width, ori_height, ori_width + 1, ori_height, ori_width, ori_height + 1 };
  int * corner_pts[4] = { wolfsburg_pts, jena_pts, stuttgart_pts, muchen_pts };

  // gen District Code
  genDeDistrictCode(corner_pts);

  // reorganize regions by district
  reorganizeRegions(charRegionArray);

  for(int i=0; i<4; ++i) {
    regions[i].bullshitsFilter();
  }

}

char PageDeutschland::judgeNPE(const double& A, const double& B, const double& C, const double& x, const double& y) {
  double res = A * x + B * y + C;

  if (res > 0) {
    return 'P';
  } else if (res < 0) {
    return 'N';
  } else {
    return 'E';
  }

}

PageDistrict::PageDistrict() {
  p_charRegionArray = NULL;
}

void PageDistrict::setRegionArray(vector<CharRegion>& charRegionArray, double g_A_base_line, double g_B_base_line,
    double g_C_base_line, double g_A_base_line_vert, double g_B_base_line_vert, double g_C_base_line_vert) {
  p_charRegionArray = &charRegionArray;
  A_base_line = g_A_base_line;
  B_base_line = g_B_base_line;
  C_base_line = g_C_base_line;
  A_base_line_vert = g_A_base_line_vert;
  B_base_line_vert = g_B_base_line_vert;
  C_base_line_vert = g_C_base_line_vert;

  bool success;
  CharRegion::getIntersectionPoint(A_base_line, B_base_line, C_base_line, A_base_line_vert, B_base_line_vert,
      C_base_line_vert, x_center, y_center, success);
}

int PageDeutschland::encode2NPE(char NPE1, char NPE2) {
  int res = 0x00;
  if (NPE1 == 'N') {
    res = 0x00;
  } else if (NPE1 == 'P') {
    res = 0x01;
  } else {
    return -1;
  }
  if (NPE2 == 'N') {
    res = res | 0x00;
  } else if (NPE2 == 'P') {
    res = res | 0x10;
  } else {
    return -1;
  }
  return res;
}

void PageDeutschland::genDeDistrictCode(int ** corner_pts) {
  int codes[4];
  for (int de_i = 0; de_i < 4; ++de_i) {
    int * corner_pt_ax = corner_pts[de_i];
    for (int cand_i = 0; cand_i < 3; ++cand_i) {
      int x = corner_pt_ax[cand_i * 2 + 0];
      int y = corner_pt_ax[cand_i * 2 + 1];
      char NPE1 = judgeNPE(A_base_line, B_base_line, C_base_line, x, y);
      char NPE2 = judgeNPE(A_base_line_vert, B_base_line_vert, C_base_line_vert, x, y);
      int code = encode2NPE(NPE1, NPE2);
      if (code != -1) {
        codes[de_i] = code;
        break;
      }
    }
  }
  code_wolfsburg = codes[0];
  code_jena = codes[1];
  code_stuttgart = codes[2];
  code_munchen = codes[3];
}

void PageDeutschland::reorganizeRegions(vector<CharRegion>& charRegionArray) {
  using namespace cv;

  Mat canvas(Size(ori_width, ori_height), CV_8UC3, CharRegion::solarized_palette["base3"]);
  line(canvas, Point(x1_base_line, y1_base_line), Point(x2_base_line, y2_base_line),
      CharRegion::solarized_palette["magenta"], 3);
  line(canvas, Point(x1_base_line_vert, y1_base_line_vert), Point(x2_base_line_vert, y2_base_line_vert),
      CharRegion::solarized_palette["magenta"], 3);

//	Mat wolfsburg_im = canvas.clone();
//	Mat jena_im = canvas.clone();
//	Mat stuttgart_im = canvas.clone();
//	Mat munchen_im = canvas.clone();
  wolfsburg_im = canvas.clone();
  jena_im = canvas.clone();
  stuttgart_im = canvas.clone();
  munchen_im = canvas.clone();

  for (int i_region = 0; i_region < charRegionArray.size(); ++i_region) {
    CharRegion& region = charRegionArray[i_region];
    for (int i_line = 0; i_line < region.const_nlines; ++i_line) {
      int x_midpoint = (region.bb_intersection_pts[i_line * 4 + 0] + region.bb_intersection_pts[i_line * 4 + 2]) / 2;
      int y_midpoint = (region.bb_intersection_pts[i_line * 4 + 1] + region.bb_intersection_pts[i_line * 4 + 3]) / 2;
      char NPE_code1 = judgeNPE(A_base_line, B_base_line, C_base_line, x_midpoint, y_midpoint);
      char NPE_code2 = judgeNPE(A_base_line_vert, B_base_line_vert, C_base_line_vert, x_midpoint, y_midpoint);
      int code = encode2NPE(NPE_code1, NPE_code2);

      // switch(code) {
      // case code_wolfsburg:
      if (code == code_wolfsburg) {
        p_wolfsburg->lines_array.push_back(make_pair(i_region, i_line));
        // break;
        line(wolfsburg_im,
            Point(region.bb_intersection_pts[i_line * 4 + 0], region.bb_intersection_pts[i_line * 4 + 1]),
            Point(region.bb_intersection_pts[i_line * 4 + 2], region.bb_intersection_pts[i_line * 4 + 3]),
            CharRegion::solarized_palette["green"], 2);
      }
      // case code_jena:
      else if (code == code_jena) {
        p_jena->lines_array.push_back(make_pair(i_region, i_line));
        // break;
        line(jena_im, Point(region.bb_intersection_pts[i_line * 4 + 0], region.bb_intersection_pts[i_line * 4 + 1]),
            Point(region.bb_intersection_pts[i_line * 4 + 2], region.bb_intersection_pts[i_line * 4 + 3]),
            CharRegion::solarized_palette["green"], 2);
      }
      // case code_stuttgart:
      else if (code == code_stuttgart) {
        p_stuttgart->lines_array.push_back(make_pair(i_region, i_line));
        // break;
        line(stuttgart_im,
            Point(region.bb_intersection_pts[i_line * 4 + 0], region.bb_intersection_pts[i_line * 4 + 1]),
            Point(region.bb_intersection_pts[i_line * 4 + 2], region.bb_intersection_pts[i_line * 4 + 3]),
            CharRegion::solarized_palette["green"], 2);
      }
      // case code_munchen:
      else if (code == code_munchen) {
        p_munchen->lines_array.push_back(make_pair(i_region, i_line));
        // break;
        line(munchen_im, Point(region.bb_intersection_pts[i_line * 4 + 0], region.bb_intersection_pts[i_line * 4 + 1]),
            Point(region.bb_intersection_pts[i_line * 4 + 2], region.bb_intersection_pts[i_line * 4 + 3]),
            CharRegion::solarized_palette["green"], 2);
      }

    }
  }

  imwrite("district_wolfsburg.png", wolfsburg_im);
  imwrite("district_jena_im.png", jena_im);
  imwrite("district_stuttgart_im.png", stuttgart_im);
  imwrite("district_munchen_im.png", munchen_im);

  for (int i = 0; i < 4; ++i) {
    regions[i].getSlopeArray(north_unit_vector);
  }
}

int g_count = 0;

void PageDistrict::getSlopeArray(double* north_unit_vector) {

  assert(p_charRegionArray != NULL);
  assert(lines_array.size() > 0);	  // !! THIS WON'T PASS IF THERE IS NO REGION IN THIS DISTRICT

  id_data_array = vector<pair<int, vector<double>>>();
  int count_pos = 0;
  int count_neg = 0;

  for (int i_line = 0; i_line < lines_array.size(); ++i_line) {
    int region_id = lines_array[i_line].first;
    int line_id = lines_array[i_line].second;
    double A = ((*p_charRegionArray)[region_id]).line_paramsABC_global[line_id * 3 + 0];
    double B = ((*p_charRegionArray)[region_id]).line_paramsABC_global[line_id * 3 + 1];
    double C = ((*p_charRegionArray)[region_id]).line_paramsABC_global[line_id * 3 + 2];

    double x, y;
    bool success;
    CharRegion::getIntersectionPoint(A, B, C, A_base_line_vert, B_base_line_vert, C_base_line_vert, x, y, success);

    // cout << g_count << "  " << x-x_center << " , " << y-y_center << endl;
    double line_dist_x = x - x_center;
    double line_dist_y = y - y_center;
    double x_vec = line_dist_x / north_unit_vector[0];
    double y_vec = line_dist_y / north_unit_vector[1];

    // if diff is too large(>0.01), won't push into id_data_array
    // cout << g_count << " : " << x_vec << " , " << y_vec << endl;
    if (abs(x_vec - y_vec) / abs(x_vec + y_vec) < 0.01 / 2. || north_unit_vector[0] == 0.
        || north_unit_vector[1] == 0.) {
      double dist_value = (x_vec + y_vec) / 2.;
      if (north_unit_vector[0] == 0.) {
        dist_value = y_vec;
      }
      if (north_unit_vector[1] == 0.) {
        dist_value = x_vec;
      }
      double angle, slope;
      CharRegion::get_slope_angle(A, B, C, angle, slope);

      id_data_array.push_back(make_pair(i_line, vector<double>( { dist_value, angle, x, y })));
//      cout << g_count << " : " << dist_value << "\t" << angle << "\t";
//      cout << x_vec << "\t" << north_unit_vector[0] << "\t" << line_dist_x << "\t;";
//      cout << y_vec << "\t" << north_unit_vector[1] << "\t" << line_dist_y << endl;

      if (dist_value < 0) {
        ++count_neg;
      } else {
        ++count_pos;
      }

    }

  }

  // reorganize id_data_array by dist_value
  if (count_neg > count_pos) {
    for (int i = 0; i < id_data_array.size(); ++i) {
      id_data_array[i].second[0] = 0 - id_data_array[i].second[0];
    }
  }
  sort(id_data_array.begin(), id_data_array.end(),
      [](const pair<int, vector<double>>& a, const pair<int, vector<double>>& b) -> bool
      {
        return a.second[0]<b.second[0];
      });

  for (int i = 0; i < id_data_array.size(); ++i) {
    cout << g_count << ":\t" << id_data_array[i].first << ",\t" << id_data_array[i].second[0];
    cout << ",\t" << id_data_array[i].second[1] / M_PI * 180. << endl;
  }

  ++g_count;

}

void PageDistrict::bullshitsFilter() {
  for(int i=0; i<id_data_array.size(); ++i) {
    if(id_data_array[i].second[0]>-200 && id_data_array[i].second[0]<1000) {
      id_data_array2.push_back(id_data_array[i]);
    }
    else {
      cout << "filtered" << endl;
    }
  }
}