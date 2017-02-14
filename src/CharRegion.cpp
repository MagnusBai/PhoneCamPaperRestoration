#include "CharRegion.h"

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

#include <mrpt/math/CVectorTemplate.h>

#include <algorithm>

#include <limits>
#include <cmath>

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::gui;
using namespace mrpt::math;
using namespace mrpt::random;
using namespace mrpt::poses;
using namespace std;

mrpt::gui::CDisplayWindow3DPtr win;

// key function 1
bool CharRegion::is_measurable() {
  int w = x_max - x_min + 1;
  int h = y_max - y_min + 1;
  bool is_ok = false;

  // maybe these rule shold be set by
  if (w > 100 || h > 100) {		// large type
    is_ok = true;
  }
  if (w > 20 && h > 20 && (float(w) / float(h) > 2.f || float(h) / float(w) > 2.f)) {	// small type
    is_ok = true;
  }

  return is_ok;
}

unordered_map<string, cv::Scalar> CharRegion::solarized_palette = {
    { string("base03"), cv::Scalar(54, 43, 0) },		// dark bg
    { string("base02"), cv::Scalar(117, 110, 88) }, { string("base01"), cv::Scalar(66, 54, 7) },
    { string("base3"), cv::Scalar(227, 246, 253) },		// light bg
    { string("yellow"), cv::Scalar(0, 137, 181) }, { string("orange"), cv::Scalar(22, 75, 203) }, {
        string("magenta"), cv::Scalar(130, 54, 211) },
    { string("violet"), cv::Scalar(196, 113, 108) }, { string("blue"), cv::Scalar(210, 139, 38) }, {
        string("cyan"), cv::Scalar(152, 161, 42) }, { string("green"), cv::Scalar(0, 153, 133) }, };

// key function 2
void CharRegion::ransac_find_lines() {
  if (count_edges < 30) {
    cout << "num of edge pixes if not enought!!!" << endl;
    return;
  }

  CVectorDouble xs, ys;

  // using char edge
//	for (int i = 0; i < count_edges; ++i) {
//		xs.push_back(char_edge[i * 2] - x_min);
//		ys.push_back(char_edge[i * 2 + 1] - y_min);
//	}

  // using contours edge
  for (int i = 0; i < count_contours; ++i) {
    xs.push_back(region_contours[i * 2] - x_min);
    ys.push_back(region_contours[i * 2 + 1] - y_min);
  }

  vector<pair<size_t, TLine2D> > detectedLines;
  const double DIST_THRESHOLD = 1.;

  ransac_detect_2D_lines(xs, ys, detectedLines, DIST_THRESHOLD, 20);

  // Display output:
  cout << "RANSAC method: ransac_detect_2D_lines" << endl;
  cout << " " << detectedLines.size() << " lines detected." << endl;

  // SAVE LINE PARAMS
  const_nlines = detectedLines.size();
  cout << "const_nlines: " << const_nlines << endl << endl;
  line_paramsABC = vector<double>(3 * const_nlines, 0.);
  line_paramsPts = vector<int>(4 * const_nlines, 0);
  line_pts = vector<int>(const_nlines, 0);

  // for(int i=0; i<const_nlines; ++i) {
  int i = 0;
  for (vector<pair<size_t, TLine2D> >::iterator p = detectedLines.begin(); p != detectedLines.end();
      ++p) {

    int count_pt = p->first;
    double A = p->second.coefs[0];
    double B = p->second.coefs[1];
    double C = p->second.coefs[2];
    int x1, y1, x2, y2;
    ABC_2points(A, B, C, x1, y1, x2, y2);
    line_pts[i] = count_pt;
    line_paramsABC[i * 3 + 0] = A;
    line_paramsABC[i * 3 + 1] = B;
    line_paramsABC[i * 3 + 2] = C;
    line_paramsPts[i * 4 + 0] = x1;
    line_paramsPts[i * 4 + 1] = y1;
    line_paramsPts[i * 4 + 2] = x2;
    line_paramsPts[i * 4 + 3] = y2;
    ++i;
  }

  // convert local points to global's
  cvtLocal2Global();

  calc_lines_info();

  calcIntersectionLine();

  vector<vector<int>> id_hierachy;
  basicHierarchicalAlg<float>(line_angles, M_PI / 60., id_hierachy);		// 180/60=3

  line_clusters_ids = id_hierachy;

  // get region_angles & region_angles_score
  for (int i = 0; i < id_hierachy.size(); ++i) {
    float angle = line_angles[id_hierachy[i][0]];
    float score = 0.f;
    if (id_hierachy[i].size() > 1) {
      score = 1.f * id_hierachy[i].size();	// angle scoring rules
    } else {
      score = 1.f;
    }
    region_angle_scores.push_back(make_pair(angle, score));
  }
}

CharRegion::CharRegion(const int arg_count_pixs, const int arg_count_contours,
    const int arg_count_edges, const int arg_x_max, const int arg_x_min, const int arg_y_max,
    const int arg_y_min) :
    count_pixs(arg_count_pixs), count_contours(arg_count_contours), count_edges(arg_count_edges), x_max(
        arg_x_max), x_min(arg_x_min), y_max(arg_y_max), y_min(arg_y_min), icap_pixs(0), icap_contours(
        0), icap_edges(0) {

  region_pixs.resize(count_pixs * 2, 0);
  region_contours.resize(count_contours * 2, 0);
  char_edge.resize(count_edges * 2, 0);

  const_nlines = 0;
  line_paramsABC = vector<double>(0);
  line_paramsPts = vector<int>(0);
  line_paramsABC_global = vector<double>(0);
  line_paramsPts_global = vector<int>(0);
  bb_intersection_pts = vector<int>(0);
  line_pts = vector<int>(0);

  line_slopes = vector<float>(0);
  line_angles = vector<float>(0);

  region_angle_scores = vector<pair<float, float>>(0);
}

void CharRegion::push_region_pix(const int arg_x, const int arg_y) {
  region_pixs[2 * icap_pixs] = arg_x;
  region_pixs[2 * icap_pixs + 1] = arg_y;
  ++icap_pixs;
}

void CharRegion::push_region_contours(const int arg_x, const int arg_y) {
  region_contours[2 * icap_contours] = arg_x;
  region_contours[2 * icap_contours + 1] = arg_y;
  ++icap_contours;
}

void CharRegion::push_edges(const int arg_x, const int arg_y) {
  char_edge[2 * icap_edges] = arg_x;
  char_edge[2 * icap_edges + 1] = arg_y;
  ++icap_edges;
}

void CharRegion::plot_char_region(const string& filename) {
  cout << y_max << " , " << y_min << endl;
  cout << x_max << " , " << x_min << endl;
  cv::Mat im(y_max - y_min + 1, x_max - x_min + 1, CV_8UC3, CharRegion::solarized_palette["base3"]);

  // draw region_pix
  for (int i = 0; i < count_pixs; ++i) {
    int x = region_pixs[i * 2] - x_min;
    int y = region_pixs[i * 2 + 1] - y_min;
    // im.at < uchar > (y, x) = 50;
    im.at<cv::Vec3b>(y, x)[0] = CharRegion::solarized_palette["violet"][0];
    im.at<cv::Vec3b>(y, x)[1] = CharRegion::solarized_palette["violet"][1];
    im.at<cv::Vec3b>(y, x)[2] = CharRegion::solarized_palette["violet"][2];
  }

  // draw edge
  for (int i = 0; i < count_edges; ++i) {
    int x = char_edge[i * 2] - x_min;
    int y = char_edge[i * 2 + 1] - y_min;
    im.at<cv::Vec3b>(y, x)[0] = CharRegion::solarized_palette["base03"][0];
    im.at<cv::Vec3b>(y, x)[1] = CharRegion::solarized_palette["base03"][1];
    im.at<cv::Vec3b>(y, x)[2] = CharRegion::solarized_palette["base03"][2];
  }

  // draw line
  for (int i = 0; i < const_nlines; ++i) {
    cv::line(im, cv::Point(line_paramsPts[i * 4 + 0], line_paramsPts[i * 4 + 1]),
        cv::Point(line_paramsPts[i * 4 + 2], line_paramsPts[i * 4 + 3]),
        CharRegion::solarized_palette["magenta"], 1);
  }

  cv::imwrite(filename.c_str(), im);
}

void CharRegion::ABC_2points(const double A, const double B, const double C, int& x1, int& y1,
    int& x2, int& y2, double X_LO, double X_UP, double Y_LO, double Y_UP) {
  double x1d, y1d, x2d, y2d;
  ABC_2points(A, B, C, x1d, y1d, x2d, y2d, X_LO, X_UP, Y_LO, Y_UP);
  y1 = int(round(y1d));
  y2 = int(round(y2d));
  x1 = int(round(x1d));
  x2 = int(round(x2d));
}

void CharRegion::ABC_2points(const double A, const double B, const double C, double& x1, double& y1,
    double& x2, double& y2, double X_LO, double X_UP, double Y_LO, double Y_UP) {
  // A x + B y + C = 0
  if (B == 0.) {		// x = xxx
    x1 = -C / A;
    x2 = x1;
    y1 = Y_LO;
    y2 = Y_UP;
  } else if (A == 0.) {		// y = yyy
    x1 = X_LO;
    x2 = X_UP;
    y1 = -C / B;
    y2 = y1;
  } else {
    double base = X_UP;	// setting larger than most image is ok
    x1 = X_LO;
    y1 = -C / B;
    x2 = base;
    y2 = -(C + A * base) / B;
  }
}

void CharRegion::twoPoints_ABC(const double x1, const double y1, const double x2, const double y2,
    double& A, double& B, double& C) {
  double move = x2 - x1;
  double rise = y2 - y1;
  if (move == 0) {
    A = 1;
    B = 0;
    C = -x1;
  } else {
    double slope = rise / move;
    double temp = y1 - slope * x1;
    A = slope;
    B = -1.;
    C = temp;
  }
}

void CharRegion::twoPoints_ABC(const int x1, const int y1, const int x2, const int y2, double& A,
    double& B, double& C) {
  double x1_d = x1;
  double y1_d = y1;
  double x2_d = x2;
  double y2_d = y2;
  CharRegion::twoPoints_ABC(x1_d, y1_d, x2_d, y2_d, A, B, C);
}

// omit the procedure of checking whether x, y  in the given line
// omit the procedure of checking (A==0||B==0)
void CharRegion::getVerticalLine(const double A, const double B, const double C,
    const int intersection_x, const int intersection_y, double& A_vert, double& B_vert,
    double& C_vert) {
  // A x + By +C = 0
  if (B == 0) {		// given line: x = xxx
    // res line: y - intersection_y = 0
    A_vert = 0.;
    B_vert = 1.;
    C_vert = -intersection_y;
  } else if (A == 0) {		// given line: y = yyy
    // res line: x - intersection_x = 0
    A_vert = 1.;
    B_vert = 0.;
    C_vert = -intersection_x;
  } else {
    double slope = -A / B;
    double slope_new = -1. / slope;
    double temp = intersection_y - slope_new * intersection_x;
    A_vert = slope_new;
    B_vert = -1.;
    C_vert = temp;
  }

}

void CharRegion::cvtLocal2Global() {
  line_paramsABC_global = line_paramsABC;
  line_paramsPts_global = line_paramsPts;
  for (int i = 0; i < const_nlines; ++i) {
    double A_p = line_paramsABC_global[i * 3 + 0];
    double B_p = line_paramsABC_global[i * 3 + 1];
    double C_p = line_paramsABC_global[i * 3 + 2] - (A_p * x_min + B_p * y_min);
    line_paramsABC_global[i * 3 + 2] = C_p;
    int x1_p, y1_p, x2_p, y2_p;
    ABC_2points(A_p, B_p, C_p, x1_p, y1_p, x2_p, y2_p);
    line_paramsPts_global[i * 4 + 0] = x1_p;
    line_paramsPts_global[i * 4 + 1] = y1_p;
    line_paramsPts_global[i * 4 + 2] = x2_p;
    line_paramsPts_global[i * 4 + 3] = y2_p;
  }
}

void CharRegion::getIntersectionPoint(const double A1, const double B1, const double C1,
    const double A2, const double B2, const double C2, int& x, int& y, bool& flag) {
  double x_d = x;
  double y_d = y;
  getIntersectionPoint(A1, B1, C1, A2, B2, C2, x_d, y_d, flag);
  x = round(x_d);
  y = round(y_d);
}

void CharRegion::getIntersectionPoint(const double A1, const double B1, const double C1,
    const double A2, const double B2, const double C2, double& x, double& y, bool& flag) {

  double slope1 = -A1 / B1;
  double slope2 = -A2 / B2;
  if (slope1 == slope2) {
    x = -numeric_limits<double>::infinity();
    y = x;
    flag = false;
    return;
  }

  y = -(C1 * A2 - C2 * A1) / (B1 * A2 - B2 * A1);
  x = -(C1 * B2 - C2 * B1) / (A1 * B2 - A2 * B1);
  flag = true;

  return;
}

// slope = rise/move = y1-y0/x1-x0
// angle is in the interval of [0, pi] radians
void CharRegion::get_slope_angle(double A, double B, double C, float& slope, float& angle) {
  slope = -float(A) / float(B);
//  angle = atan(1.f / slope);	// atan(x), in the
//  // interval [-pi/2,+pi/2] radians.
//  if (angle < 0) {
//    angle += M_PI;
//  }
  angle = get_angle_by_slope(slope);
}

void CharRegion::get_slope_angle(double A, double B, double C, double& slope, double& angle) {
  float slope_f, angle_f;
  get_slope_angle(A, B, C, slope_f, angle_f);
  slope = slope_f;
  angle = angle_f;
}

double CharRegion::get_angle_by_slope(const double slope) {
  double angle = atan(1. / slope);    // atan(x), in the
  // interval [-pi/2,+pi/2] radians.
  if (angle < 0.) {
    angle += M_PI;
  }
  return angle;
}

float CharRegion::get_angle_by_slope(const float slope) {
  double slope_d = slope;
  double angle_d = get_angle_by_slope(slope_d);
  return static_cast<float>(angle_d);
}

//float CharRegion::get_angle_diff(float a1, float a2) {
//	float d1 = abs(a1-a2);
//	float d2 = M_PI-d1;
//	return min(d1, d2);
//}
template<typename T>
T CharRegion::get_angle_diff(T a1, T a2) {
  T d1 = abs(a1 - a2);
  T d2 = M_PI - d1;
  return min(d1, d2);
}

void CharRegion::calc_lines_info() {
  line_slopes.resize(const_nlines);
  line_angles.resize(const_nlines);

// calc each lines angle and slope
  for (int i = 0; i < const_nlines; ++i) {
    double A = line_paramsABC_global[i * 3 + 0];
    double B = line_paramsABC_global[i * 3 + 1];
    double C = line_paramsABC_global[i * 3 + 2];
    cout << "(" << A << " , " << B << " , " << C << ")" << endl;
    float slope, angle;
    get_slope_angle(A, B, C, slope, angle);
    line_slopes[i] = slope;
    line_angles[i] = angle;
  }

// cluster lines in an very closed angle
}

// ----------
// DESIGN PATTERN: BETTER PASS DIST_FUN() AS ARGUMENTS
template<typename T>
void basicHierarchicalAlg(const vector<T>& data, const T diff_thresh,
    vector<vector<int>>& id_hierachy) {
// temp printing
  cout << endl << endl << "data:~~ thresh:" << diff_thresh << endl;
  for (int i = 0; i < data.size(); ++i) {
    cout << "  (" << i << ")" << data[i] << "  ";
  }
  cout << endl << endl;
// temp printing

  int count_data = data.size();
  vector<vector<T>> diff_mat = vector<vector<T>>(count_data, vector<T>(count_data, T(0)));
  for (int h = 0; h < count_data; ++h) {
    for (int w = 0; w < count_data; ++w) {
      //   EXTRA COMPUTATION , symmetrical matrix may save 1/2 time
      diff_mat[h][w] = CharRegion::get_angle_diff<T>(data[h], data[w]);
    }
  }

// init remaining_data_dict
  unordered_map<int, int> remaining_data_dict;
  for (int i = 0; i < count_data; ++i) {
    remaining_data_dict[i] = 0;
  }
// setting
  unordered_map<int, int> todo_data_dict;
  vector<int> group_ids;

  while (remaining_data_dict.size() != 0) {
    int start_id = remaining_data_dict.begin()->first;								// extract first elem
    todo_data_dict[start_id] = 0;
    remaining_data_dict.erase(start_id);

    group_ids.clear();

    while (todo_data_dict.size() != 0) {
      int cur_id = todo_data_dict.begin()->first;
      todo_data_dict.erase(cur_id);
      group_ids.push_back(cur_id);
      // for(int id=cur_id+1; id<count_data; ++id) {
      for (int id = 0; id < count_data; ++id) {
        if (id == cur_id) {
          continue;
        }
        // I think (cur_id, id) can void equivalent pair of <a, b>, <b, a>
        if (diff_mat[cur_id][id] <= diff_thresh) {
          if (remaining_data_dict.find(id) != remaining_data_dict.end()) {
            // in the remaining_data_dict
            todo_data_dict[id] = 0;
            remaining_data_dict.erase(id);
          }
        }		// is an elem that can be grouped together
      }		// id: [cur_id+1: count_data]
    }		// todo_data_dict.size!=0

    sort(group_ids.begin(), group_ids.end());

    id_hierachy.push_back(group_ids);

  }		// remaining_data_dict.size()!=0

  sort(id_hierachy.begin(), id_hierachy.end(),
      [](const vector<int>& a, const vector<int>& b) -> bool
      {
        return a.size() > b.size();
      });

// temp printing
  for (int i_cluster = 0; i_cluster < id_hierachy.size(); ++i_cluster) {
    vector<int>& ids = id_hierachy[i_cluster];
    cout << " [  ";
    for (int i_elems = 0; i_elems < ids.size(); ++i_elems) {
      cout << "(" << ids[i_elems] << ")" << data[ids[i_elems]] << "  ";
    }
    cout << "] \n";
  }
// temp printing
}

void CharRegion::calcIntersectionLine() {
  bb_intersection_pts = vector<int>(const_nlines * 4, 0);
  for (int i = 0; i < const_nlines; ++i) {
    float A = line_paramsABC_global[i * 3 + 0];
    float B = line_paramsABC_global[i * 3 + 1];
    float C = line_paramsABC_global[i * 3 + 2];
    vector<float> v_params;		// x1, y1, x2, y2;

    // aux var
    float x_, y_;
    // #1
    x_ = x_min;
    y_ = -(A * x_ + C) / B;
    if (y_ >= y_min && y_ <= y_max) {
      v_params.push_back(int(round(x_)));
      v_params.push_back(int(round(y_)));
    }
    // #2
    x_ = x_max;
    y_ = -(A * x_ + C) / B;
    if (y_ >= y_min && y_ <= y_max) {
      v_params.push_back(int(round(x_)));
      v_params.push_back(int(round(y_)));
    }
    // #3
    y_ = y_min;
    x_ = -(B * y_ + C) / A;
    if (x_ >= x_min && x_ <= x_max) {
      v_params.push_back(int(round(x_)));
      v_params.push_back(int(round(y_)));
    }
    // #4
    y_ = y_max;
    x_ = -(B * y_ + C) / A;
    if (x_ >= x_min && x_ <= x_max) {
      v_params.push_back(int(round(x_)));
      v_params.push_back(int(round(y_)));
    }

    if (v_params.size() == 4) {
      bb_intersection_pts[i * 4 + 0] = v_params[0];
      bb_intersection_pts[i * 4 + 1] = v_params[1];
      bb_intersection_pts[i * 4 + 2] = v_params[2];
      bb_intersection_pts[i * 4 + 3] = v_params[3];
    } else {
      bb_intersection_pts[i * 4 + 0] = -1;
      bb_intersection_pts[i * 4 + 1] = -1;
      bb_intersection_pts[i * 4 + 2] = -1;
      bb_intersection_pts[i * 4 + 3] = -1;
    }
  }
}

//template<typename T>
//void clac1DHistogram(const vector<T>& data, const int count_hist, const pair<T, T>& range,
//		vector<int>& hist, vector<vector<int>>& hist_ids) {
//	hist = vector<int>(count_hist, 0);
//	hist_ids = vector<vector<int>>(count_hist, vector<int>(0));
//	T cell_range=(range.secon+static_cast<T>(0.001f)-range.first)/static_cast<T>(count_hist);
//	for( int i=0; i<data.size(); ++i ) {
//		int id = static_cast<int>( floor(data[i]/cell_range) );
//		++hist[id];
//		hist_ids[id].push_back(i);
//	}
//}

//template <> void clac1DHistogram<float>(const vector<float>& data, const int count_hist,
//		const pair<float, float>& range, vector<int>& hist, vector<vector<int>>& hist_ids);
