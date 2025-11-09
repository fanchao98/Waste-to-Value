#pragma once
#ifndef SAMPLE_ON_BALL_
#define SAMPLE_ON_BALL_

#include <vector>
#include <cmath>

#include "helpers.h"
#include <opencv2/imgproc/imgproc.hpp>
#define PI 3.1415926
class SAMPLE_ON_BALL
{
public:
    void OrientationSamplePoints();
    void OrientationSamplePoints_2();
    std::vector<cv::Point3d> sample_points;
};


#endif 