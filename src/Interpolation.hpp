#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;

// 三线性插值
double trilinearInterpolation(const Vector3d& point, const std::vector<Vector3d>& vertices, const std::vector<double>& values)
{
    double x = point.x();
    double y = point.y();
    double z = point.z();

    double c = 0.0;

    for (int i = 0; i < 8; i++) {
        double weight = 1.0;
        for (int j = 0; j < 3; j++) {
            double v = (i >> j) & 1 ? 1 - point[j] : point[j];
            weight *= v;
        }
        c += weight * values[i];
    }

    return c;
}

