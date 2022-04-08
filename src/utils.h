
#pragma once

#include "misc_v14.h"

//------------------------------------------------
bool collision_test_point_ellipse(double x, double y, double f1x, double f1y,
                                  double f2x, double f2y, double e);
//------------------------------------------------
bool collision_test_abline_line(double m, double k,
                                double l1x, double l1y, double l2x, double l2y,
                                double &hx, double &hy);
//------------------------------------------------
bool collision_test_abline_ellipse(double m, double k,
                                   double f1x, double f1y, double f2x, double f2y, double e,
                                   double &h1x, double &h1y, double &h2x, double &h2y);

//------------------------------------------------
bool collision_test_line_ellipse(double l1x, double l1y, double l2x, double l2y,
                                 double f1x, double f1y, double f2x, double f2y, double e);

//------------------------------------------------
bool collision_test_circle_rect(double rx, double ry, double w, double h, double theta,
                                double cx, double cy, double r);

//------------------------------------------------
bool collision_test_hex_ellipse(double hx, double hy, double w, double f1x,
                                double f1y, double f2x, double f2y, double e);

