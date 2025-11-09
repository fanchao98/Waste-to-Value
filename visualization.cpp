#include"visualization.h"

void Visual::insert_Line(General_Mesh& mesh, std::vector<Vec3> points, float radius)
{
	if (points.size() < 2)
		return;

	for (int i = 0; i < points.size() - 1; i++)
	{

		Vec3 v0 = points[i];
		Vec3 v1 = points[i + 1];

		Vec3 t_dir = v1 - v0;

		Vec3 scale(radius, radius, t_dir.Length());
		Vec3 trans = (v0 + v1) / 2;

		vector<Vec3> c_vertex = mesh.meshScale(scale, mesh.cylinder_vertex);
		c_vertex = mesh.meshRotate(t_dir, c_vertex);
		c_vertex = mesh.meshTrans(trans, c_vertex);
		mesh.insert(c_vertex, mesh.cylinder_faces);
	}
}

void Visual::creat_red_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points)
{
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;

	igl::readOBJ("ball.obj", V_2, F_2);
	ofstream all_balls(file_name + "_unaccessivle_points.obj");
	for (int i = 0; i < vis_points.size(); i++) {
		for (int j = 0; j < V_2.rows(); j++)
			all_balls << "v " << V_2(j, 0) + vis_points[i](0, 0) << " " << V_2(j, 1) + vis_points[i](1, 0) << " " << V_2(j, 2) + vis_points[i](2, 0) << " 0.9" << " 0.05" << " 0.05" << endl;
		for (int j = 0; j < F_2.rows(); j++)
			all_balls << "f " << F_2(j, 0) + i * V_2.rows() + 1 << " " << F_2(j, 1) + i * V_2.rows() + 1 << " " << F_2(j, 2) + i * V_2.rows() + 1 << endl;
	}
	all_balls.close();
}

void Visual::creat_red_ball_2(string file_name, std::vector<Eigen::MatrixXd> vis_points)
{
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;

	igl::readOBJ("ball.obj", V_2, F_2);
	ofstream all_balls(file_name + "_unaccessivle_points.obj");
	for (int i = 0; i < vis_points.size(); i++) {
		for (int j = 0; j < V_2.rows(); j++)
			all_balls << "v " << V_2(j, 0) + vis_points[i](0, 0) << " " << V_2(j, 1) + vis_points[i](1, 0) << " " << V_2(j, 2) + vis_points[i](2, 0) << " 0.7" << " 0.15" << " 0.15" << endl;
		for (int j = 0; j < F_2.rows(); j++)
			all_balls << "f " << F_2(j, 0) + i * V_2.rows() + 1 << " " << F_2(j, 1) + i * V_2.rows() + 1 << " " << F_2(j, 2) + i * V_2.rows() + 1 << endl;
	}
	all_balls.close();
}

void Visual::creat_green_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points, vector<double> color_map)
{
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;
	
	igl::readOBJ("ball.obj", V_2, F_2);
	ofstream all_balls(file_name + "_covering_points.obj");
	for (int i = 0; i < vis_points.size(); i++) {
		for (int j = 0; j < V_2.rows(); j++)
			all_balls << "v " << V_2(j, 0) + vis_points[i](0, 0) << " " << V_2(j, 1) + vis_points[i](1, 0) << " " << V_2(j, 2) + vis_points[i](2, 0) << " " << "0.05" << " " << color_map[i] * 1.0 << " " << "0.05" << endl;
		for (int j = 0; j < F_2.rows(); j++)
			all_balls << "f " << F_2(j, 0) + i * V_2.rows() + 1 << " " << F_2(j, 1) + i * V_2.rows() + 1 << " " << F_2(j, 2) + i * V_2.rows() + 1 << endl;
	}
	all_balls.close();
}

void Visual::generateModelForRendering(vector<Eigen::Vector3d> lines, string file_name, std::vector<Eigen::MatrixXd> vis_points)
{
	float radius = 0.5;
	General_Mesh mesh1;
	mesh1.r = 0.7;
	mesh1.g = 0.15;
	mesh1.b = 0.15;
	std::vector<Vec3> points;

	Vec3 the_point;

	Vec3 the_point_2;
	for (int i = 0; i < vis_points.size(); i++) {
		the_point_2.m_x = vis_points[i](0, 0);
		the_point_2.m_y = vis_points[i](1, 0);
		the_point_2.m_z = vis_points[i](2, 0);

		the_point.m_x = the_point_2.m_x + lines[i].x() * 2;
		the_point.m_y = the_point_2.m_y + lines[i].y() * 2;
		the_point.m_z = the_point_2.m_z + lines[i].z() * 2;
		points.push_back(the_point_2);
		points.push_back(the_point);
		insert_Line(mesh1, points, radius);
		points.clear();
	}
	mesh1.genResultMesh(("temp_vis\\" + file_name + ".obj").c_str());
	points.clear();
}