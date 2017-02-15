/*
 * PagePartition.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: magnus
 */

#include <algorithm>

#include "PagePartition.h"
#include <list>
#include <iterator>

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

  for (int i = 0; i < 4; ++i) {
    regions[i].bullshitsFilter();
    regions[i].oddDistanceFilter(3, regions[i].id_data_array2);
  }

  for (int i = 0; i < p_wolfsburg->id_data_array2.size(); ++i) {
    int id = p_wolfsburg->id_data_array2[i].first;
    int region_id = p_wolfsburg->lines_array[id].first;
    int line_id = p_wolfsburg->lines_array[id].second;
    int x1 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 0];
    int y1 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 1];
    int x2 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 2];
    int y2 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 3];

    line(wolfsburg_im, Point(x1, y1), Point(x2, y2), CharRegion::solarized_palette["violet"], 2);
  }
  imwrite("district_wolfsburg_im2.png", wolfsburg_im);

  for (int i = 0; i < p_wolfsburg->id_data_array3.size(); ++i) {
    int id = p_wolfsburg->id_data_array3[i].first;
    int region_id = p_wolfsburg->lines_array[id].first;
    int line_id = p_wolfsburg->lines_array[id].second;
    int x1 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 0];
    int y1 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 1];
    int x2 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 2];
    int y2 = (p_wolfsburg->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 3];

    line(wolfsburg_im, Point(x1, y1), Point(x2, y2), CharRegion::solarized_palette["magenta"], 2);
  }
  imwrite("district_wolfsburg_im3.png", wolfsburg_im);

  for (int i = 0; i < p_jena->id_data_array3.size(); ++i) {
    int id = p_jena->id_data_array3[i].first;
    int region_id = p_jena->lines_array[id].first;
    int line_id = p_jena->lines_array[id].second;
    int x1 = (p_jena->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 0];
    int y1 = (p_jena->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 1];
    int x2 = (p_jena->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 2];
    int y2 = (p_jena->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 3];

    line(jena_im, Point(x1, y1), Point(x2, y2), CharRegion::solarized_palette["magenta"], 2);
  }
  imwrite("district_jena_im3.png", jena_im);

  for (int i = 0; i < p_munchen->id_data_array3.size(); ++i) {
    int id = p_munchen->id_data_array3[i].first;
    int region_id = p_munchen->lines_array[id].first;
    int line_id = p_munchen->lines_array[id].second;
    int x1 = (p_munchen->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 0];
    int y1 = (p_munchen->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 1];
    int x2 = (p_munchen->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 2];
    int y2 = (p_munchen->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 3];

    line(munchen_im, Point(x1, y1), Point(x2, y2), CharRegion::solarized_palette["magenta"], 2);
  }
  imwrite("district_munchen_im3.png", munchen_im);

  for (int i = 0; i < p_stuttgart->id_data_array3.size(); ++i) {
    int id = p_stuttgart->id_data_array3[i].first;
    int region_id = p_stuttgart->lines_array[id].first;
    int line_id = p_stuttgart->lines_array[id].second;
    int x1 = (p_stuttgart->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 0];
    int y1 = (p_stuttgart->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 1];
    int x2 = (p_stuttgart->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 2];
    int y2 = (p_stuttgart->p_charRegionArray)->at(region_id).bb_intersection_pts[line_id * 4 + 3];

    line(stuttgart_im, Point(x1, y1), Point(x2, y2), CharRegion::solarized_palette["magenta"], 2);
  }
  imwrite("district_stuttgart_im3.png", stuttgart_im);

  for (int i = 0; i < 4; ++i) {
    regions[i].tellIfSheered();
  }

  cout << x1_base_line_vert << " , " << y1_base_line_vert << " , " << x2_base_line_vert << " , " << y2_base_line_vert << endl;
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

  circle(canvas, Point(x1_base_line, y1_base_line), 10, CharRegion::solarized_palette["violet"], 3);
  circle(canvas, Point(x2_base_line, y2_base_line), 10, CharRegion::solarized_palette["violet"], 3);
  circle(canvas, Point(x1_base_line_vert, y1_base_line_vert), 10, CharRegion::solarized_palette["violet"], 3);
  circle(canvas, Point(x2_base_line_vert, y2_base_line_vert), 10, CharRegion::solarized_palette["violet"], 3);

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
    cout << ",\t" << id_data_array[i].second[1] << endl;
  }

  ++g_count;

}

void PageDistrict::bullshitsFilter() {
  for (int i = 0; i < id_data_array.size(); ++i) {
    if (id_data_array[i].second[0] > -200 && id_data_array[i].second[0] < 1000) {
      id_data_array2.push_back(id_data_array[i]);
    } else {
      cout << "filtered" << endl;
    }
  }
}

void PageDistrict::oddDistanceFilter(int scan_range, vector<pair<int, vector<double>>>& source_signal,
double odd_angle_thresh,
int count_odd_thresh) {
  id_data_array3.clear();
  list<pair<int, vector<double>>> s_signal_list(source_signal.begin(), source_signal.end());

  for(list<pair<int, vector<double>>>::iterator iter=s_signal_list.begin(); iter!=s_signal_list.end(); ++iter) {
    vector<double> pre_vec(scan_range, 0.);
    vector<double> pro_vec(scan_range, 0.);
    bool WOULD_CALC = true;

    list<pair<int, vector<double>>>::iterator it_probe;
    // traversal backwards
    it_probe = iter;
    for(int i=0; i<scan_range; ++i) {
      if( it_probe==s_signal_list.begin() ) {
        WOULD_CALC = false;
        break;
      }
      --it_probe;
      pre_vec[i] = (*it_probe).second[1];
    }
    // traversal forwards
    it_probe = iter;
    ++it_probe;
    for(int i=0; i<scan_range; ++i) {
      if( it_probe==s_signal_list.end() ) {
        WOULD_CALC = false;
        break;
      }
      pro_vec[i] = (*it_probe).second[1];
      ++it_probe;
    }
    if(WOULD_CALC) {
      double error_thresh = odd_angle_thresh/180.*M_PI;
      int count_odd_values = 0;
      double this_angle = (*iter).second[1];

      for(int i=0; i<3; ++i) {
        double angle1 = pre_vec[i];
        double angle2 = pro_vec[i];
        if( CharRegion::get_angle_diff<double>(angle1, this_angle)>=error_thresh ) {
          ++count_odd_values;
        }
        if( CharRegion::get_angle_diff<double>(angle2, this_angle)>=error_thresh ) {
          ++count_odd_values;
        }
      }

      for(int i=0; i<3; ++i) {
        cout << pre_vec[i] << "\t";
      }
      cout << "(" << this_angle << ")\t";
      for(int i=0; i<3; ++i) {
        cout << pro_vec[i] << "\t";
      }
      cout << endl;

      if(count_odd_values>=count_odd_thresh) {
        cout << "ODD VALUE: ";
      }
      else {
        id_data_array3.push_back(*iter);
      }

    }

  }
}

bool PageDistrict::tellIfSheered(double mean_sheered_angle) {
  float sum_angles = 0.f;

  for (int i = 0; i < id_data_array3.size(); ++i) {
    int id = id_data_array3[i].first;
    int region_id = lines_array[id].first;
    int line_id = lines_array[id].second;
    float line_angle = (p_charRegionArray)->at(region_id).line_angles[line_id];   // MAYBE EXISTS SOME ERROR!!!!!!!!
    sum_angles += line_angle;
    cout << line_angle << "\t";
  }
  double angle_aver = sum_angles / id_data_array3.size();
  double base_angle, base_slope;
  CharRegion::get_slope_angle(A_base_line, B_base_line, C_base_line, base_slope, base_angle);

  cout << "\nBASE: " << base_angle << "  AVERAGE: " << angle_aver << " THRESH: " << mean_sheered_angle << endl;

  return abs(angle_aver - base_angle) >= mean_sheered_angle;
}

// according to the limit theory of advanced mathematics, when theta is small, lim( tan(theta)-theta ) is incline to 0
//0   angle:  0.0 tan():  0.0
//1   angle:  0.017 tan():  0.017
//2   angle:  0.034 tan():  0.034
//3   angle:  0.052 tan():  0.052
//4   angle:  0.069 tan():  0.069
//5   angle:  0.087 tan():  0.087
//6   angle:  0.104 tan():  0.105
//7   angle:  0.122 tan():  0.122
//8   angle:  0.139 tan():  0.140
//9   angle:  0.157 tan():  0.158
//10  angle:  0.174 tan():  0.176
//11  angle:  0.191 tan():  0.194
//12  angle:  0.209 tan():  0.212
//13  angle:  0.226 tan():  0.230
//14  angle:  0.244 tan():  0.249
//15  angle:  0.261 tan():  0.267
//16  angle:  0.279 tan():  0.286
//17  angle:  0.296 tan():  0.305
//18  angle:  0.314 tan():  0.324
//19  angle:  0.331 tan():  0.344
//20  angle:  0.349 tan():  0.363
//21  angle:  0.366 tan():  0.383
//22  angle:  0.383 tan():  0.404
//23  angle:  0.401 tan():  0.424
//24  angle:  0.418 tan():  0.445
//25  angle:  0.436 tan():  0.466
//26  angle:  0.453 tan():  0.487
//27  angle:  0.471 tan():  0.509
//28  angle:  0.488 tan():  0.531
//29  angle:  0.506 tan():  0.554
//30  angle:  0.523 tan():  0.577
