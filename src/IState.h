#pragma once

//#include"ofxMSAmcts.h"
#include <Eigen/Dense>
#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>
#include <random>
using namespace std;

//#include "ofMain.h"

//using namespace msa::mcts;


//--------------------------------------------------------------
//--------------------------------------------------------------
// 行为结构
struct Action {
    Point_3 new_normal;
    double offset;
};

#define kNumActions    30  //原来50  //30

//--------------------------------------------------------------
//--------------------------------------------------------------
// 状态类
class State {
public:

    //--------------------------------------------------------------
    // 必须有的方法（接口）

    State() {
        reset();
    }

    // 是否为终止状态（到达末尾）
    bool is_terminal() const {
        if (index_of_current_plane == total_num_of_planes)  //index_of_current_plane == total_num_of_planes - 1
            return true;
        else
            return false;
    }

    bool is_terminal_2() const {
        if (index_of_current_plane + 1 == total_num_of_planes)  //index_of_current_plane == total_num_of_planes - 1
            return true;
        else
            return false;
    }

    // 将动作应用于状态
    void apply_action(const Action& action) {
        normal_of_planes.push_back(action.new_normal);
        origin_of_planes.push_back(Point_3(barycenter_of_planes[index_of_current_plane].x() + action.new_normal.x() * action.offset,
            barycenter_of_planes[index_of_current_plane].y() + action.new_normal.y() * action.offset, barycenter_of_planes[index_of_current_plane].z() + action.new_normal.z() * action.offset));
        index_of_current_plane++;
    }

    // 返回即将做决策的代理的代理ID（从零开始）
    int agent_id() const {
        return 0;
    }

    // 返回从此状态可能的动作
    void get_actions(std::vector<Action>& actions, vector<vector<Point_3>> all_vertices_in_components, int num_first_plane_components, vector<Point_3> new_current_points)  {
        actions.resize(kNumActions);
        int step = kNumActions / 6;
        for (int i = 0; i < kNumActions; i++) {
            //依次插值五个方向 
            int index = floor(i / step);
            double xx, yy, zz;
            Vector3 v;
            if (i >= 25)
            {
                xx = selected_planes_normal[index_of_current_plane].x();
                yy = selected_planes_normal[index_of_current_plane].y();
                zz = selected_planes_normal[index_of_current_plane].z();
            }
            else if (i % 6 == 0) {
                xx = -10.0f + static_cast<float>(rand()) / RAND_MAX * 10.0;
                yy = -10.0f + static_cast<float>(rand()) / RAND_MAX * 10.0;
                zz = 0.0f + static_cast<float>(rand()) / RAND_MAX * 10.0;
            }
            v = Vector3(xx, yy, zz);
            v.Normalized();
            actions[i].new_normal = Point_3(v.m_x,v.m_y,v.m_z);
            

            //if (selected_planes_normal.size() != 0)  //将期望方向统一改为(0,0,1)
            //    actions[i].new_normal = Point_3(selected_planes_normal[index_of_current_plane].x() + (0 - selected_planes_normal[index_of_current_plane].x()) * (4 - index) / 4,
            //        selected_planes_normal[index_of_current_plane].y() + (0 - selected_planes_normal[index_of_current_plane].y()) * (4 - index) / 4,
            //        selected_planes_normal[index_of_current_plane].z() + (1 - selected_planes_normal[index_of_current_plane].z()) * (4 - index) / 4);

            int offset_ratio = i % step;
            //actions[i].new_normal = Point_3(0, 0, 1);  //暂时测试
            get_random_action(actions[i], all_vertices_in_components, num_first_plane_components, new_current_points, offset_ratio, step);
        }
    }

    double truncatedGaussianProbability(double x, double mean, double stdDev) {
        // 计算标准化横坐标值
        double z = (x - mean) / stdDev;

        // 截断边界
        double lowerBound = (0 - mean) / stdDev;
        double upperBound = (1 - mean) / stdDev;

        // 计算截断后的概率值
        double probability = (std::erf(upperBound / sqrt(2)) - std::erf(lowerBound / sqrt(2))) / 2;
        double truncatedCoefficient = 1.0 / probability;

        // 根据截断系数调整概率密度
        double exponent = -0.5 * pow(z, 2);
        double coefficient = 1.0 / (stdDev * sqrt(2 * 3.1415926));
        return coefficient * exp(exponent) * truncatedCoefficient;
    }

    double inverseErfcApproximation(double p) {
        // 近似逆误差函数的算法实现
        double t = sqrt(-2 * log(p / 2.0));
        double c0 = 2.515517;
        double c1 = 0.802853;
        double c2 = 0.010328;
        double d1 = 1.432788;
        double d2 = 0.189269;
        double d3 = 0.001308;

        return t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1);
    }

    double generateRandomX(double mean, double stdDev) {
        // 创建随机数生成器
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // 生成服从均匀分布的随机数
        double u = dis(gen);

        // 使用逆变换采样获取横坐标值
        double x = mean + stdDev * sqrt(2) * inverseErfcApproximation(2 * u);

        // 限制横坐标取值范围为0~1
        if (x < 0) {
            x = 0.01;
        }
        if (x > 1) {
            x = 0.99;
        }

        return x;
    }

    // 获取随机动作，如果没有找到动作则返回false
    bool get_random_action(Action& action, vector<vector<Point_3>> all_vertices_in_components, int num_first_plane_components, vector<Point_3> new_current_points, int offset_ratio, int step)  {
        /*std::vector<Point_3> sample_normals = OrientationSamplePoints_normal_2();
        int random_index = rand() % sample_normals.size(); 
        action.new_normal = sample_normals[random_index];*/
        
        //action.new_normal = selected_planes_normal[index_of_current_plane];  //暂时先不改变normal，取消启发式方法

        double max_offset_distance;
        Point_3 origin_point(barycenter_of_planes[index_of_current_plane].x(),barycenter_of_planes[index_of_current_plane].y(), barycenter_of_planes[index_of_current_plane].z());
        if (index_of_additive_interface[index_of_current_plane].size() == 0) {  //求最大距离
            double max_distance = -100000000;
            if (index_of_current_plane < num_first_plane_components) {
                for (int i = 0; i < all_vertices_in_components[index_of_current_plane].size(); i++) {
                    double distance = Calculate_distance_plane_and_point(action.new_normal, origin_point, all_vertices_in_components[index_of_current_plane][i]);
                    if (distance > max_distance) {
                        max_distance = distance;
                    }
                }
            }
            else {
                for (int i = 0; i < new_current_points.size(); i++) {
                    double distance = Calculate_distance_plane_and_point(action.new_normal, origin_point, new_current_points[i]);
                    if (distance > max_distance) {
                        max_distance = distance;
                    }
                }
            }
            max_offset_distance = max_distance;
        }
        else {  //求最小距离
            double min_distance = 100000000;
			for (int i = 0; i < index_of_additive_interface[index_of_current_plane].size(); i++) {
                for (int j = 0; j < additive_interface_points[index_of_additive_interface[index_of_current_plane][i]].size(); j++) {
                    double distance = Calculate_distance_plane_and_point(action.new_normal, origin_point, additive_interface_points[index_of_additive_interface[index_of_current_plane][i]][j]);
                    if (distance < min_distance) {
                        min_distance = distance;
                    }
                }
			}
            max_offset_distance = min_distance;
        }
        //cout << max_offset_distance << endl;

        double mean = (1 - mean_scalar_value_of_cutting_planes[index_of_current_plane]) * 1;
        double stdDev = 0.4;
        double random_offset_factor = generateRandomX(mean, stdDev);    //目前效果不好，起反作用   //用(当前场值-最小场值)/(最大-最小场值)来控制排序？
        double probability = truncatedGaussianProbability(random_offset_factor, mean, stdDev);  //概率计算得似乎有问题

        //float random_offset = random_offset_factor * max_offset_distance;    //目前启发式效果不好
        /*if (index_of_current_plane == 1)
            max_offset_distance = 0.1;*/
        float random_offset = 0.0f + static_cast<float>(rand()) / RAND_MAX * max_offset_distance;  //朴素方法
        //float random_offset =  max_offset_distance-0.1;    //伪随机数
        //float random_offset = 0.0f + offset_ratio * (max_offset_distance-0.1) / (step-1);
        action.offset = random_offset;
        return true;
    }

    // 评估此状态并返回奖励的向量（对于每个代理）
    const vector<float> evaluate(bool flag_satisfy_constraints,double new_volume,bool is_max_volume,double last_max_volume) {
        vector<float> rewards(1);
        V_current = new_volume;

        if(last_max_volume != 0)
            rewards[0] = 15 * (V_current - V_offset) / (V_ori - V_offset) * pow((V_current / last_max_volume),3);  //前面的常数可能需要根据不同的模型进行调整  //a=30
        else
            rewards[0] = 15 * (V_current - V_offset) / (V_ori - V_offset);
        //rewards[0] = (V_current - V_offset) / (V_ori - V_offset);
        //rewards[0] = 10;
        
        if(flag_satisfy_constraints == false) 
			rewards[0] = 0; 
        //else  //暂时测试用
        //    rewards[0] = 0.00000001;

        return rewards;
    }
        
    void reset() {
        normal_of_planes.clear();
        origin_of_planes.clear();
        V_current = 0;
        index_of_current_plane = 0;
    }

    std::vector<Point_3> OrientationSamplePoints_normal_2()
    {
        std::vector<Point_3> sample_points_normal;
        int number_points = 50;
        for (int i = 1; i <= number_points; i++) {
            double phi = acos(-1.0 + (2.0 * i - 1.0) / number_points);
            double theta = sqrt(number_points * Pi) * phi;
            double x = cos(theta) * sin(phi);
            double y = sin(theta) * sin(phi);
            double z = cos(phi);
            if (z > 0)
                sample_points_normal.push_back(Point_3(x, y, z));
        }
        return sample_points_normal;
    }
    //--------------------------------------------------------------
        
    double Calculate_distance_plane_and_point(Point_3 normal, Point_3 origin, Point_3 p)
    {
        double a = normal.x();
        double b = normal.y();
        double c = normal.z();
        double d = -a * origin.x() - b * origin.y() - c * origin.z();
        double min_distance = 100000000;
        int index_closest_point = -1;
        double distance = (a * p.x() + b * p.y() + c * p.z() + d) / std::sqrt(a * a + b * b + c * c);
        return distance;
    }

    vector<Point_3> normal_of_planes;
    vector<Point_3> origin_of_planes;
    //double current_volume;
    int index_of_current_plane;
    vector<Eigen::Vector3d> expected_normals_of_cutting_planes;
    vector<double> mean_scalar_value_of_cutting_planes;
    std::vector<std::vector<Point_3>> additive_interface_points;
    vector<vector<int>> index_of_additive_interface;


    vector<Point_3> selected_planes_normal;
    vector<Point_3> barycenter_of_planes;
    int total_num_of_planes;
    int total_num_of_complete_planes;
    double V_ori;
    double V_offset;
    double V_current;


   /* std::vector<Eigen::MatrixXd> green_points;
    std::vector<double> color_map;*/
};
