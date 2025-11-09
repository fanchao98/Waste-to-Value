#include "sample_on_ball.h"


void SAMPLE_ON_BALL::OrientationSamplePoints()
{
    sample_points.clear();
    int number_points = num_ori_sample;
    for (int i = 1; i <= number_points; i++) {
        double phi = acos(-1.0 + (2.0 * i - 1.0) / number_points);
        double theta = sqrt(number_points * PI) * phi;
        double x = cos(theta) * sin(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(phi);
        //sample_points.push_back(cv::Point3d(x, y, z));  
        if (z >= -0.2)   //暂时也限制减材方向  //z >= -0.3
            sample_points.push_back(cv::Point3d(x, y, z));
    }

}

void SAMPLE_ON_BALL::OrientationSamplePoints_2()
{
    sample_points.clear();
    int number_points = num_ori_sample / 2;  //增加竖直方向  /2
    sample_points.push_back(cv::Point3d(0, 0, 1));
    //sample_points.push_back(cv::Point3d(0.22, 0.22, 0.8));
    /*sample_points.push_back(cv::Point3d(0.22, -0.22, 0.8));
    sample_points.push_back(cv::Point3d(-0.22, 0.22, 0.8));
    sample_points.push_back(cv::Point3d(-0.22, -0.22, 0.8));*/
    for (int i = 1; i <= number_points; i++) { 
        double phi = acos(-1.0 + (2.0 * i - 1.0) / number_points);
        double theta = sqrt(number_points * PI) * phi;
        double x = cos(theta) * sin(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(phi);
        if(z >= 0)  //z >= -0.3
            sample_points.push_back(cv::Point3d(x, y, z));
    }
}
