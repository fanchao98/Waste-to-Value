#pragma once
#include<string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include<igl/readOBJ.h>
#include<igl/readSTL.h>
#include"GeneralMesh.h"
using namespace std;

class Visual
{
public:
	void insert_Line(General_Mesh& mesh, std::vector<Vec3> points, float radius);
	void creat_red_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points);
	void creat_red_ball_2(string file_name, std::vector<Eigen::MatrixXd> vis_points);
	void creat_green_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points, vector<double> color_map);
	void generateModelForRendering(vector<Eigen::Vector3d> lines, string file_name, std::vector<Eigen::MatrixXd> vis_points);
};