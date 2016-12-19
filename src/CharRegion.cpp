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

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::gui;
using namespace mrpt::math;
using namespace mrpt::random;
using namespace mrpt::poses;
using namespace std;

mrpt::gui::CDisplayWindow3DPtr win;

void CharRegion::log() {
	cout << "fuckyou" << endl;
}

void CharRegion::ransac_find_lines() {
	if(count_edges<30) {
		cout << "num of edge pixes if not enought!!!" << endl;
		return;
	}

	CVectorDouble xs, ys;

	for (int i = 0; i < count_edges; ++i) {
		xs.push_back(char_edge[i * 2]);
		ys.push_back(char_edge[i * 2 + 1]);
	}

	vector < pair<size_t, TLine2D> > detectedLines;
	const double DIST_THRESHOLD = 1.;

	CTicTac tictac;

	ransac_detect_2D_lines(xs, ys, detectedLines, DIST_THRESHOLD, 20);

	// Display output:
	cout << "RANSAC method: ransac_detect_2D_lines" << endl;
	cout << " Computation time: " << tictac.Tac() * 1000.0 << " ms" << endl;
	cout << " " << detectedLines.size() << " lines detected." << endl;

	for(vector < pair<size_t, TLine2D> >::iterator it=detectedLines.begin();
			it<detectedLines.end(); ++it) {
		cout << "\t\t" << "!!!!!" << it->first << "\t";
		cout << it->second.coefs[0] << "\t" << it->second.coefs[1] << "\t" << it->second.coefs[2];
		cout <<endl;
	}
}
