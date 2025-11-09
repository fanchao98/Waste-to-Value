#include "PPS.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

void PPS::Calculate_additive_blocks(string meshname, vector<Plane>final_cutting_planes, vector<vector<Point_3>>points_in_cutting_planes)
{
	ifstream load_final_cutting_planes("MCTS_temp\\" + mesh_target + "\\final_cutting_planes.txt");
	double a, b, c, d, e, f;
	while (load_final_cutting_planes >> a >> b >> c >> d >> e >> f) {
		Plane temp_plane(Point_3(a, b, c), Point_3(d, e, f));
		final_cutting_planes.push_back(temp_plane);
		vector<Point_3> temp_vec;
		points_in_cutting_planes.push_back(temp_vec);
		load_final_cutting_planes >> a >> b >> c;
		if (a == 0 && b == 0) {
			for (int i = 0; i < c; i++) {
				double t1, t2, t3;
				load_final_cutting_planes >> t1 >> t2 >> t3;
				points_in_cutting_planes[points_in_cutting_planes.size() - 1].push_back(Point_3(t1, t2, t3));
			}
		}
		else
			points_in_cutting_planes[points_in_cutting_planes.size() - 1].push_back(Point_3(a, b, c));
	}

	string path = "postprocess\\" + mesh_target;
	_mkdir(path.c_str());

	Slicer_2 slicer;
	Eigen::MatrixXd V_ff;
	Eigen::MatrixXi F_ff;
	igl::readOFF(meshname, V_ff, F_ff);
	igl::writeOBJ("postprocess\\" + mesh_target + "\\temp_target_mesh.obj", V_ff, F_ff);

	slicer.clear();
	slicer.load("postprocess\\" + mesh_target + "\\temp_target_mesh.obj");
	vector<int> index_exist_planes;
	index_exist_planes.clear();

	for (int i = 0; i < final_cutting_planes.size(); i++) {
		Slicer_2 removed_triangles_from_slicer;
		vector<vector<int>> cutting_plane_polygons;
		int ori_num_vertex = slicer.positions.size();
		slicer.normal[0] = final_cutting_planes[i].normal.x();
		slicer.normal[1] = final_cutting_planes[i].normal.y();
		slicer.normal[2] = final_cutting_planes[i].normal.z();
		slicer.origin[0] = final_cutting_planes[i].origin.x();
		slicer.origin[1] = final_cutting_planes[i].origin.y();
		slicer.origin[2] = final_cutting_planes[i].origin.z();
		slicer.cut();
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vector_after(final_cutting_planes[i].normal.x(), final_cutting_planes[i].normal.y(), final_cutting_planes[i].normal.z());
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();

		Slicer_2 temp_slicer = slicer;
		for (int j = 0; j < slicer.positions.size(); j++) {
			Eigen::MatrixXd temp_V;
			temp_V.resize(3, 1);
			temp_V(0, 0) = slicer.positions[j][0];
			temp_V(1, 0) = slicer.positions[j][1];
			temp_V(2, 0) = slicer.positions[j][2];
			temp_V = rotMatrix.inverse() * temp_V;
			temp_slicer.positions[j][0] = temp_V(0, 0);
			temp_slicer.positions[j][1] = temp_V(1, 0);
			temp_slicer.positions[j][2] = temp_V(2, 0);
		}
		//temp_slicer.save("postprocess\\" + mesh_target + "\\aaaaa.obj");
		Eigen::MatrixXd temp_cutting_point;
		temp_cutting_point.resize(3, 1);
		temp_cutting_point(0, 0) = final_cutting_planes[i].origin.x();
		temp_cutting_point(1, 0) = final_cutting_planes[i].origin.y();
		temp_cutting_point(2, 0) = final_cutting_planes[i].origin.z();
		temp_cutting_point = rotMatrix.inverse() * temp_cutting_point;

		//save triangles above the cutting plane
		removed_triangles_from_slicer.clear();
		removed_triangles_from_slicer.positions = slicer.positions;
		vector<int> temp_index;
		vector<Point_3> points_in_cutting_planes_boundary;
		for (int j = 0; j < slicer.triangles.size(); j++) {
			for (int k = 0; k < 3; k++) {
				if (abs(temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0)) < 0.00001 || abs(temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0)) < 0.00001
					|| abs(temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0)) < 0.00001) {
					rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vectorBefore).toRotationMatrix();
					Eigen::MatrixXd temp_V_2;
					temp_V_2.resize(3, 1);
					temp_V_2(0, 0) = temp_slicer.positions[slicer.triangles[j][k]][0];
					temp_V_2(1, 0) = temp_slicer.positions[slicer.triangles[j][k]][1];
					temp_V_2(2, 0) = temp_slicer.positions[slicer.triangles[j][k]][2];
					temp_V_2 = rotMatrix.inverse() * temp_V_2;
					points_in_cutting_planes_boundary.push_back(Point_3(temp_V_2(0, 0), temp_V_2(1, 0), temp_V_2(2, 0)));
				}

			}
			for (int k = 0; k < 3; k++) {
				if (temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0) > 0.0000001) {
					removed_triangles_from_slicer.triangles.push_back(slicer.triangles[j]);
					temp_index.push_back(j);
					//slicer.triangles.erase(slicer.triangles.begin() + j);
					//j--;
					break;
				}
			}
		}
		vector<double> min_distance(points_in_cutting_planes[i].size());
		vector<int> min_index(points_in_cutting_planes[i].size());
		for (int j = 0; j < points_in_cutting_planes[i].size(); j++) {
			min_distance[j] = 100000;
			min_index[j] = -1;
		}
		for (int j = 0; j < points_in_cutting_planes_boundary.size(); j++) {
			for (int k = 0; k < points_in_cutting_planes[i].size(); k++) {
				double temp_distance = CGAL::squared_distance(points_in_cutting_planes_boundary[j], points_in_cutting_planes[i][k]);
				if (temp_distance < min_distance[k]) {
					min_distance[k] = temp_distance;
					min_index[k] = j;
				}
			}
		}

		removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\temp_removed_mesh.obj");

		Mesh mesh;
		CGAL::IO::read_OBJ("postprocess\\" + mesh_target + "\\temp_removed_mesh.obj", mesh);
		std::vector<std::vector<int>> connected_components_triangles;
		connected_components_triangles.clear();
		vector<bool> flag_face_been_search(mesh.number_of_faces());
		for (int j = 0; j < mesh.number_of_faces(); j++)
			flag_face_been_search[j] = false;
		typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
		const double bound = std::cos(0.75 * CGAL_PI);
		std::vector<face_descriptor> cc;
		vector<double> min_distance_2(points_in_cutting_planes[i].size());
		vector<int> min_index_2(points_in_cutting_planes[i].size());
		for (int j = 0; j < points_in_cutting_planes[i].size(); j++) {
			min_distance_2[j] = 100000;
			min_index_2[j] = -1;
		}
		for (face_descriptor fd : faces(mesh))
		{
			cc.clear();
			for (Mesh::Vertex_index vd : mesh.vertices_around_face(mesh.halfedge(fd))) {
				Point_3 temp_point = mesh.point(vd);
				for (int k = 0; k < points_in_cutting_planes[i].size(); k++) {
					double temp_distance = CGAL::squared_distance(temp_point, points_in_cutting_planes_boundary[min_index[k]]);
					if (temp_distance < min_distance_2[k]) {
						min_distance_2[k] = temp_distance;
						min_index_2[k] = fd.idx();
					}
				}
			}
			if (flag_face_been_search[fd.idx()] == true)
				continue;
			PMP::connected_component(fd,
				mesh,
				std::back_inserter(cc));
			std::vector<int> temp_points;
			connected_components_triangles.push_back(temp_points);
			for (int j = 0; j < cc.size(); j++) {
				flag_face_been_search[cc[j].idx()] = true;
				connected_components_triangles[connected_components_triangles.size() - 1].push_back(cc[j].idx());
			}
		}


		vector<int> index_removed_component(points_in_cutting_planes[i].size());
		for (int j = 0; j < points_in_cutting_planes[i].size(); j++) {
			index_removed_component[j] = -1;
		}
		for (int m = 0; m < points_in_cutting_planes[i].size(); m++) {
			bool flag_break = false;
			for (int j = 0; j < connected_components_triangles.size(); j++) {
				for (int k = 0; k < connected_components_triangles[j].size(); k++) {
					//cout << j << k << m<<endl;
					if (connected_components_triangles[j][k] == min_index_2[m]) {  //是需要切除的连通分量
						flag_break = true;
						index_removed_component[m] = j;
						break;
					}
				}
				if (flag_break)
					break;
			}
		}
		Slicer_2 real_removed_triangles_from_slicer, temp_slicer2;
		real_removed_triangles_from_slicer.positions = slicer.positions;
		for (int m = 0; m < points_in_cutting_planes[i].size(); m++) {
			temp_slicer2 = slicer;
			for (int t = 0; t < connected_components_triangles[index_removed_component[m]].size(); t++) {
				real_removed_triangles_from_slicer.triangles.push_back(temp_slicer2.triangles[temp_index[connected_components_triangles[index_removed_component[m]][t]]]);
				temp_slicer2.triangles.erase(temp_slicer2.triangles.begin() + temp_index[connected_components_triangles[index_removed_component[m]][t]]);
				for (int j = t + 1; j < connected_components_triangles[index_removed_component[m]].size(); j++)
					if (temp_index[connected_components_triangles[index_removed_component[m]][j]] > temp_index[connected_components_triangles[index_removed_component[m]][t]])
						temp_index[connected_components_triangles[index_removed_component[m]][j]]--;
			}
		}
		real_removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\real_removed_mesh-" + to_string(i) + ".obj");
		vector<bool> jud_triangles_exist(slicer.triangles.size());
		for (int j = 0; j < slicer.triangles.size(); j++) {
			jud_triangles_exist[j] = true;
			for (int k = 0; k < real_removed_triangles_from_slicer.triangles.size(); k++) {
				if (real_removed_triangles_from_slicer.triangles[k][0] == slicer.triangles[j][0] && real_removed_triangles_from_slicer.triangles[k][1] == slicer.triangles[j][1] && real_removed_triangles_from_slicer.triangles[k][2] == slicer.triangles[j][2]) {
					jud_triangles_exist[j] = false;
					break;
				}
			}
		}
		temp_slicer2 = slicer;
		temp_slicer2.triangles.clear();
		for (int j = 0; j < slicer.triangles.size(); j++) {
			if (jud_triangles_exist[j])
				temp_slicer2.triangles.push_back(slicer.triangles[j]);
		}
		slicer = temp_slicer2;

		//填充洞
		Add_triangles_of_cutting_planes(slicer);
		slicer.save("postprocess\\" + mesh_target + "\\temp_target_mesh.obj");
	}

	//对每个removed mesh进行切割以满足自支撑约束
	Segmentation_for_each_removed_mesh(final_cutting_planes);

	//建图，重新计算切片层加工顺序
	Sort_slice_layers();

}

void PPS::Segmentation_for_each_removed_mesh(vector<Plane>final_cutting_planes)
{
	//目前切片不足一个步长会报错
	float step_slice = 2;    //切片步长
	double threshold_of_support;  //自支撑约束阈值  //0.1  //0.2
	double threshold_of_support_current = 1.0;  //当前层自支撑约束阈值  //0.4

	all_slice_points.clear();
	all_slice_normal.clear();
	for (int t = 0; t < final_cutting_planes.size(); t++) {
		int cont_block = 0;
		Slicer_2 slicer;
		slicer.load("postprocess\\" + mesh_target + "\\real_removed_mesh-" + to_string(t) + ".obj");
		Point_3 selected_normal;
		Point_3 selected_origin;
		//double z_cutting = -1;
		double z_ori_cutting = -1;
		bool jud_break = false;
		while (true) {
			cout << endl;
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(final_cutting_planes[t].normal.x(), final_cutting_planes[t].normal.y(), final_cutting_planes[t].normal.z());
			Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int j = 0; j < slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = slicer.positions[j][0];
				temp_V(1, 0) = slicer.positions[j][1];
				temp_V(2, 0) = slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				slicer.positions[j][0] = temp_V(0, 0);
				slicer.positions[j][1] = temp_V(1, 0);
				slicer.positions[j][2] = temp_V(2, 0);
			}
			Eigen::MatrixXd temp_V_2(3, 1);
			temp_V_2(0, 0) = final_cutting_planes[t].origin.x();
			temp_V_2(1, 0) = final_cutting_planes[t].origin.y();
			temp_V_2(2, 0) = final_cutting_planes[t].origin.z();
			temp_V_2 = rotMatrix.inverse() * temp_V_2;
			z_ori_cutting = temp_V_2(2, 0);
			//slicer.save("postprocess\\" + mesh_target + "\\ttttt.obj");
			OrientationSamplePoints();
			double max_volume = -100000;
			for (int s = 0; s < sample_points_sphere.size(); s++) {
				//selected_normal = sample_points_sphere[s];
				if (s == 0)
					threshold_of_support = 0.2;  //0.2
				else
					threshold_of_support = 0.1;  //0.1
				bool jud_no_cutting = false;
				Slicer_2 temp_slicer = slicer;
				vectorBefore = Eigen::Vector3d(0, 0, 1);
				vector_after = Eigen::Vector3d(sample_points_sphere[s].x(), sample_points_sphere[s].y(), sample_points_sphere[s].z());
				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
				for (int j = 0; j < slicer.positions.size(); j++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = slicer.positions[j][0];
					temp_V(1, 0) = slicer.positions[j][1];
					temp_V(2, 0) = slicer.positions[j][2];
					temp_V = rotMatrix.inverse() * temp_V;
					temp_slicer.positions[j][0] = temp_V(0, 0);
					temp_slicer.positions[j][1] = temp_V(1, 0);
					temp_slicer.positions[j][2] = temp_V(2, 0);
				}
				double max_z = -100000, min_z = 100000;
				for (int j = 0; j < temp_slicer.triangles.size(); j++) {
					for (int k = 0; k < 3; k++) {
						if (temp_slicer.positions[temp_slicer.triangles[j][k]][2] > max_z)
							max_z = temp_slicer.positions[temp_slicer.triangles[j][k]][2];
						else if (temp_slicer.positions[temp_slicer.triangles[j][k]][2] < min_z)
							min_z = temp_slicer.positions[temp_slicer.triangles[j][k]][2];
					}
				}

				Slicer_2 temp_slicer_2;
				Slicer_2 removed_triangles_from_slicer;
				Slicer_2 removed_triangles_from_slicer_2;
				Slicer_2 removed_triangles_from_current_layer;
				if (max_z - step_slice <= min_z) {
					jud_no_cutting = true;
					if (sample_points_sphere[s].x() == 0 && sample_points_sphere[s].y() == 0 && sample_points_sphere[s].z() == 1) {
						selected_origin = Point_3(temp_V_2(0, 0), temp_V_2(1, 0), z_ori_cutting);
						selected_normal = sample_points_sphere[s];
						temp_slicer_2 = temp_slicer;
						temp_slicer_2.save("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
						break;
					}
				}
				bool flag_break = false;
				for (float slice_z = max_z - step_slice; slice_z > min_z; slice_z -= step_slice) {
					temp_slicer_2 = temp_slicer;
					temp_slicer_2.normal = { 0,0,1 };
					temp_slicer_2.origin = { 0,0,slice_z };
					temp_slicer_2.cut();
					if (slice_z != max_z - step_slice) {
						temp_slicer_2.origin = { 0,0,slice_z + step_slice };
						temp_slicer_2.cut();
						temp_slicer_2.origin = { 0,0,slice_z };
					}

					removed_triangles_from_slicer.positions = temp_slicer_2.positions;
					removed_triangles_from_slicer.triangles.clear();
					removed_triangles_from_slicer.origin = temp_slicer_2.origin;
					removed_triangles_from_slicer.normal = { 0,0,1 };
					removed_triangles_from_current_layer.positions = temp_slicer_2.positions;
					removed_triangles_from_current_layer.triangles.clear();
					removed_triangles_from_current_layer.origin = temp_slicer_2.origin;
					removed_triangles_from_current_layer.normal = { 0,0,1 };
					for (int j = 0; j < temp_slicer_2.triangles.size(); j++) {
						for (int k = 0; k < 3; k++) {
							if (temp_slicer_2.positions[temp_slicer_2.triangles[j][k]][2] - slice_z > 0.0000001 && temp_slicer_2.positions[temp_slicer_2.triangles[j][k]][2] - slice_z < step_slice + 0.0000001) {
								removed_triangles_from_current_layer.triangles.push_back(temp_slicer_2.triangles[j]);
							}
							if (temp_slicer_2.positions[temp_slicer_2.triangles[j][k]][2] - slice_z > 0.0000001) {
								removed_triangles_from_slicer.triangles.push_back(temp_slicer_2.triangles[j]);
								temp_slicer_2.triangles.erase(temp_slicer_2.triangles.begin() + j);
								j--;
								break;
							}
						}
					}

					bool jud_flat_area = false;
					double sum_area = 0, sum_area_2 = 0;;
					double area_need_support = Calculate_area_need_support(removed_triangles_from_slicer, sum_area);
					removed_triangles_from_slicer_2 = removed_triangles_from_slicer;
					double area_need_support_current = Calculate_area_need_support_2(removed_triangles_from_current_layer, sum_area_2, jud_flat_area);

					if (area_need_support / sum_area > threshold_of_support || area_need_support_current / sum_area_2 > threshold_of_support_current || (jud_flat_area && (sample_points_sphere[s].x() != 0 || sample_points_sphere[s].y() != 0 || sample_points_sphere[s].z() != 1))) {
						//if (area_need_support / sum_area > threshold_of_support || area_need_support_current / sum_area_2 > threshold_of_support_current ) {

								//removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\temp_slice.obj");
						removed_triangles_from_slicer.load("postprocess\\" + mesh_target + "\\temp_slice.obj");
						removed_triangles_from_slicer_2 = removed_triangles_from_slicer;
						temp_slicer_2.load("postprocess\\" + mesh_target + "\\temp_slicer_remain.obj");
						flag_break = true;
						break;
					}
					//removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\temp_slice.obj");
					bool jud_break_3;
					jud_break_3 = Check_cutting_bottom(removed_triangles_from_slicer);
					if (jud_break_3 == false) {
						removed_triangles_from_slicer.load("postprocess\\" + mesh_target + "\\temp_slice.obj");
						removed_triangles_from_slicer_2 = removed_triangles_from_slicer;
						temp_slicer_2.load("postprocess\\" + mesh_target + "\\temp_slicer_remain.obj");
						flag_break = true;
						break;
					}

					removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\temp_slice.obj");
					temp_slicer_2.save("postprocess\\" + mesh_target + "\\temp_slicer_remain.obj");
					//z_cutting = slice_z;
				}
				Mesh temp_mesh;
				CGAL::IO::read_OBJ("postprocess\\" + mesh_target + "\\temp_slice.obj", temp_mesh);
				double volume = CGAL::Polygon_mesh_processing::volume(temp_mesh);
				if (volume > max_volume && !jud_no_cutting) {
					cout << volume << endl;
					max_volume = volume;
					vectorBefore = Eigen::Vector3d(sample_points_sphere[s].x(), sample_points_sphere[s].y(), sample_points_sphere[s].z());
					vector_after = Eigen::Vector3d(0, 0, 1);
					rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
					for (int j = 0; j < removed_triangles_from_slicer_2.positions.size(); j++) {
						Eigen::MatrixXd temp_V;
						temp_V.resize(3, 1);
						temp_V(0, 0) = removed_triangles_from_slicer_2.positions[j][0];
						temp_V(1, 0) = removed_triangles_from_slicer_2.positions[j][1];
						temp_V(2, 0) = removed_triangles_from_slicer_2.positions[j][2];
						temp_V = rotMatrix.inverse() * temp_V;
						removed_triangles_from_slicer_2.positions[j][0] = temp_V(0, 0);
						removed_triangles_from_slicer_2.positions[j][1] = temp_V(1, 0);
						removed_triangles_from_slicer_2.positions[j][2] = temp_V(2, 0);
					}
					double max_z_temp_slicer_2 = -100000;
					int index_max_z_temp_slicer_2 = -1;
					int index_max_z_temp_slicer_2_2 = -1;
					for (int j = 0; j < temp_slicer_2.triangles.size(); j++)
						for (int k = 0; k < 3; k++) {
							if (temp_slicer_2.positions[temp_slicer_2.triangles[j][k]][2] > max_z_temp_slicer_2) {
								max_z_temp_slicer_2 = temp_slicer_2.positions[temp_slicer_2.triangles[j][k]][2];
								index_max_z_temp_slicer_2 = j;
								index_max_z_temp_slicer_2_2 = k;
							}
						}
					for (int j = 0; j < temp_slicer_2.positions.size(); j++) {
						Eigen::MatrixXd temp_V;
						temp_V.resize(3, 1);
						temp_V(0, 0) = temp_slicer_2.positions[j][0];
						temp_V(1, 0) = temp_slicer_2.positions[j][1];
						temp_V(2, 0) = temp_slicer_2.positions[j][2];
						temp_V = rotMatrix.inverse() * temp_V;
						temp_slicer_2.positions[j][0] = temp_V(0, 0);
						temp_slicer_2.positions[j][1] = temp_V(1, 0);
						temp_slicer_2.positions[j][2] = temp_V(2, 0);
					}
					Eigen::MatrixXd temp_V_2;
					temp_V_2.resize(3, 1);
					temp_V_2(0, 0) = temp_slicer_2.positions[temp_slicer_2.triangles[index_max_z_temp_slicer_2][index_max_z_temp_slicer_2_2]][0];
					temp_V_2(1, 0) = temp_slicer_2.positions[temp_slicer_2.triangles[index_max_z_temp_slicer_2][index_max_z_temp_slicer_2_2]][1];
					temp_V_2(2, 0) = temp_slicer_2.positions[temp_slicer_2.triangles[index_max_z_temp_slicer_2][index_max_z_temp_slicer_2_2]][2];
					//temp_V_2 = rotMatrix.inverse() * temp_V_2;
					selected_origin = Point_3(temp_V_2(0, 0), temp_V_2(1, 0), temp_V_2(2, 0));
					removed_triangles_from_slicer_2.save("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
					temp_slicer_2.save("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
					selected_normal = sample_points_sphere[s];
					if (!flag_break && s == 0)
						break;
				}
			}
			if (selected_normal.x() == 0 && selected_normal.y() == 0 && selected_normal.z() == 1 && abs(selected_origin.z() - z_ori_cutting) <= step_slice) {
				Slicer_2 slicer_a, slicer_b;
				slicer_a.load("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				slicer_b.load("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				slicer_a.positions = slicer_b.positions;
				for (int i = 0; i < slicer_b.triangles.size(); i++) {
					slicer_a.triangles.push_back(slicer_b.triangles[i]);
					slicer_b.triangles.erase(slicer_b.triangles.begin() + i);
					i--;
				}
				slicer_b.save("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				//旋转回去
				Eigen::Vector3d vectorBefore(final_cutting_planes[t].normal.x(), final_cutting_planes[t].normal.y(), final_cutting_planes[t].normal.z());
				Eigen::Vector3d vector_after(0, 0, 1);
				Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
				for (int j = 0; j < slicer_a.positions.size(); j++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = slicer_a.positions[j][0];
					temp_V(1, 0) = slicer_a.positions[j][1];
					temp_V(2, 0) = slicer_a.positions[j][2];
					temp_V = rotMatrix.inverse() * temp_V;
					slicer_a.positions[j][0] = temp_V(0, 0);
					slicer_a.positions[j][1] = temp_V(1, 0);
					slicer_a.positions[j][2] = temp_V(2, 0);
				}
				slicer_a.save("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				slicer.normal = { 0,0,1 };
				jud_break = true;
			}
			if (!jud_break) {
				slicer.load("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				slicer.normal = { -selected_normal.x(),-selected_normal.y(),-selected_normal.z() };
				slicer.origin = { selected_origin.x(),selected_origin.y(),selected_origin.z() };
				Add_triangles_of_cutting_planes_2(slicer);
				slicer.normal = { selected_normal.x(),selected_normal.y(),selected_normal.z() };
				//旋转回去
				Eigen::Vector3d vectorBefore(final_cutting_planes[t].normal.x(), final_cutting_planes[t].normal.y(), final_cutting_planes[t].normal.z());
				Eigen::Vector3d vector_after(0, 0, 1);
				//Eigen::Vector3d vector_after(selected_normal.x(), selected_normal.y(), selected_normal.z());
				Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
				for (int j = 0; j < slicer.positions.size(); j++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = slicer.positions[j][0];
					temp_V(1, 0) = slicer.positions[j][1];
					temp_V(2, 0) = slicer.positions[j][2];
					temp_V = rotMatrix.inverse() * temp_V;
					slicer.positions[j][0] = temp_V(0, 0);
					slicer.positions[j][1] = temp_V(1, 0);
					slicer.positions[j][2] = temp_V(2, 0);
				}
				slicer.save("postprocess\\" + mesh_target + "\\remain_additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");

				Slicer_2 slicer_a;
				slicer_a.load("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
				for (int j = 0; j < slicer_a.positions.size(); j++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = slicer_a.positions[j][0];
					temp_V(1, 0) = slicer_a.positions[j][1];
					temp_V(2, 0) = slicer_a.positions[j][2];
					temp_V = rotMatrix.inverse() * temp_V;
					slicer_a.positions[j][0] = temp_V(0, 0);
					slicer_a.positions[j][1] = temp_V(1, 0);
					slicer_a.positions[j][2] = temp_V(2, 0);
				}
				slicer_a.save("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj");
			}

			////对每个子区域沿各自方向逐层切片得到切片层
			Eigen::MatrixXd V_ff;
			Eigen::MatrixXi F_ff;
			igl::readOBJ("postprocess\\" + mesh_target + "\\additive_block-" + to_string(t) + "-" + to_string(cont_block) + ".obj", V_ff, F_ff);

			//先换算下
			Eigen::Vector3d vectorBefore2 = Eigen::Vector3d(final_cutting_planes[t].normal.x(), final_cutting_planes[t].normal.y(), final_cutting_planes[t].normal.z());
			Eigen::Vector3d vector_after2 = Eigen::Vector3d(0, 0, 1);
			Eigen::Matrix3d rotMatrix2 = Eigen::Quaterniond::FromTwoVectors(vectorBefore2, vector_after2).toRotationMatrix();
			Eigen::MatrixXd temp_normal;
			temp_normal.resize(3, 1);
			temp_normal(0, 0) = slicer.normal[0];
			temp_normal(1, 0) = slicer.normal[1];
			temp_normal(2, 0) = slicer.normal[2];
			temp_normal = rotMatrix2.inverse() * temp_normal;
			slicer.normal[0] = temp_normal(0, 0);
			slicer.normal[1] = temp_normal(1, 0);
			slicer.normal[2] = temp_normal(2, 0);

			//按各自方向旋转
			vectorBefore = Eigen::Vector3d(0, 0, 1);
			vector_after = Eigen::Vector3d(slicer.normal[0], slicer.normal[1], slicer.normal[2]);
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int j = 0; j < V_ff.rows(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = V_ff(j, 0);
				temp_V(1, 0) = V_ff(j, 1);
				temp_V(2, 0) = V_ff(j, 2);
				temp_V = rotMatrix.inverse() * temp_V;
				V_ff(j, 0) = temp_V(0, 0);
				V_ff(j, 1) = temp_V(1, 0);
				V_ff(j, 2) = temp_V(2, 0);
			}
			igl::writeSTL("postprocess\\" + mesh_target + "\\temp_for_slice.stl", V_ff, F_ff);
			Katana::Instance().config.loadConfig("config.ini");
			Katana::Instance().stl.loadStl(("postprocess\\" + mesh_target + "\\temp_for_slice.stl").c_str());
			Katana::Instance().slicer.buildLayers();
			Katana::Instance().slicer.buildSegments();
			all_slice_points.push_back(vector<vector<Vertex>>());
			all_slice_normal.push_back(vector<Vertex>());
			all_index_slice.push_back(pair<int, int>(t, cont_block));
			Katana::Instance().gcode.write(all_slice_points[all_slice_points.size() - 1]);
			//旋转回来
			vectorBefore = Eigen::Vector3d(slicer.normal[0], slicer.normal[1], slicer.normal[2]);
			vector_after = Eigen::Vector3d(0, 0, 1);
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int j = 0; j < all_slice_points[all_slice_points.size() - 1].size(); j++) {
				all_slice_normal[all_slice_normal.size() - 1].push_back(Vertex(slicer.normal[0], slicer.normal[1], slicer.normal[2]));

				for (int k = 0; k < all_slice_points[all_slice_points.size() - 1][j].size(); k++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = all_slice_points[all_slice_points.size() - 1][j][k].x;
					temp_V(1, 0) = all_slice_points[all_slice_points.size() - 1][j][k].y;
					temp_V(2, 0) = all_slice_points[all_slice_points.size() - 1][j][k].z;
					temp_V = rotMatrix.inverse() * temp_V;
					all_slice_points[all_slice_points.size() - 1][j][k].x = temp_V(0, 0);
					all_slice_points[all_slice_points.size() - 1][j][k].y = temp_V(1, 0);
					all_slice_points[all_slice_points.size() - 1][j][k].z = temp_V(2, 0);
				}
			}
			if (jud_break)
				break;
			cont_block++;
		}
	}
}

void PPS::Sort_slice_layers()
{
	//修改尺寸，考虑误差
	Nozzle.upper_surface_r = 2.5; //2.5
	Nozzle.lowwer_surface_r = 0.2;  //0.2
	Nozzle.nozzle__H_total = 10; //10
	Nozzle.nozzle_H_half = 5;

	vector<vector<My_node>> all_nodes(all_slice_normal.size());
	for (int i = 0; i < all_slice_normal.size(); i++)
		all_nodes[i].resize(all_slice_normal[i].size());

	ofstream ttemp_file("postprocess\\" + mesh_target + "\\aaa.obj");
	for (int i = 0; i < all_slice_points.size(); i++) {
		for (int j = 0; j < all_slice_points[i].size(); j++) {
			for (int k = 0; k < all_slice_points[i][j].size(); k++)
				ttemp_file << "v " << all_slice_points[i][j][k].x << " " << all_slice_points[i][j][k].y << " " << all_slice_points[i][j][k].z << endl;
		}
	}

	//碰撞检测，建立边
	for (int sub = all_slice_points.size() - 1; sub >= 0; sub--) {
		for (int contour = 0; contour < all_slice_points[sub].size(); contour++) {
			for (int last_sub = all_slice_points.size() - 1; last_sub >= 0; last_sub--) {
				if (last_sub == sub)
					continue;
				for (int last_contour = 0; last_contour < all_slice_points[last_sub].size(); last_contour++) {
					Eigen::Vector3d vectorBefore(0, 0, 1);
					Eigen::Vector3d vector_after(all_slice_normal[sub][contour].x, all_slice_normal[sub][contour].y, all_slice_normal[sub][contour].z);
					Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
					bool jud_collision = false;
					for (int p = 0; p < all_slice_points[sub][contour].size(); p++) {
						Eigen::MatrixXd temp_V;
						temp_V.resize(3, 1);
						temp_V(0, 0) = all_slice_points[sub][contour][p].x;
						temp_V(1, 0) = all_slice_points[sub][contour][p].y;
						temp_V(2, 0) = all_slice_points[sub][contour][p].z;
						temp_V = rotMatrix.inverse() * temp_V;
						for (int last_p = 0; last_p < all_slice_points[last_sub][last_contour].size(); last_p += 2) {
							Eigen::MatrixXd temp_V_2;
							temp_V_2.resize(3, 1);
							temp_V_2(0, 0) = all_slice_points[last_sub][last_contour][last_p].x;
							temp_V_2(1, 0) = all_slice_points[last_sub][last_contour][last_p].y;
							temp_V_2(2, 0) = all_slice_points[last_sub][last_contour][last_p].z;
							temp_V_2 = rotMatrix.inverse() * temp_V_2;
							double circle_r;
							if (temp_V_2(2, 0) < temp_V(2, 0))
								continue;
							else if (temp_V_2(2, 0) - temp_V(2, 0) < Nozzle.nozzle_H_half)
								circle_r = Nozzle.lowwer_surface_r + (temp_V_2(2, 0) - temp_V(2, 0)) * (Nozzle.upper_surface_r - Nozzle.lowwer_surface_r) / Nozzle.nozzle_H_half;
							else if (temp_V_2(2, 0) - temp_V(2, 0) < Nozzle.nozzle__H_total)
								circle_r = Nozzle.upper_surface_r;
							else {
								jud_collision = true;
								break;
							}
							if (pow(all_slice_points[last_sub][last_contour][last_p].x - all_slice_points[sub][contour][p].x, 2) + pow(all_slice_points[last_sub][last_contour][last_p].y - all_slice_points[sub][contour][p].y, 2)
								- pow(circle_r, 2) < 0) {
								jud_collision = true;
								break;
							}
						}
						if (jud_collision)
							break;
					}
					if (jud_collision == true) {   //添加边
						all_nodes[sub][contour].edges.push_back(make_pair(last_sub, last_contour));
						//cout << sub << " " << contour << " " << last_sub << " " << last_contour << endl;
						all_nodes[last_sub][last_contour].indegree++;
						//break;
					}
				}
				/*if (jud_collision == true)
					break;*/
			}
		}
	}

	//在图中搜索，得到切片层加工顺序
	new_slicer_order.clear();
	vector<vector<bool>> flag_node_been_search(all_nodes.size());
	for (int i = 0; i < all_nodes.size(); i++) {
		flag_node_been_search[i].resize(all_nodes[i].size());
		for (int j = 0; j < all_nodes[i].size(); j++)
			flag_node_been_search[i][j] = false;
	}
	bool jud_break = false;
	while (true) {
		bool jud_all_search = true;
		for (int sub = all_nodes.size() - 1; sub >= 0; sub--) {
			for (int contour = 0; contour < all_nodes[sub].size(); contour++) {
				if (all_nodes[sub][contour].indegree > 0) {
					jud_all_search = false;
					break;
				}
				if (all_nodes[sub][contour].indegree == 0 && flag_node_been_search[sub][contour] == false) {
					new_slicer_order.push_back(make_pair(sub, contour));
					flag_node_been_search[sub][contour] = true;
					for (int i = 0; i < all_nodes[sub][contour].edges.size(); i++) {
						all_nodes[all_nodes[sub][contour].edges[i].first][all_nodes[sub][contour].edges[i].second].indegree--;
					}
					jud_all_search = false;
				}
				if (sub == 0 && contour == all_nodes[sub].size() - 1 && jud_all_search == true) {
					jud_break = true;
					break;
				}
			}
			if (jud_break)
				break;
		}
		if (jud_break)
			break;
	}

	ofstream out("postprocess\\" + mesh_target + "\\test_vis.obj");
	int next_j = -1, next_k = -1;
	int cont_blocks = 0;
	vector<bool> jud_slicer(all_slice_points.size());
	for (int i = 0; i < all_slice_points.size(); i++)
		jud_slicer[i] = false;
	for (int i = 0; i < new_slicer_order.size(); i++)
		jud_slicer[new_slicer_order[i].first] = true;
	for (int i = 0; i < all_slice_points.size(); i++) {
		if (jud_slicer[i] == false) {
			int j = all_index_slice[i].first;
			int k = all_index_slice[i].second;
			Slicer_2 slicer;
			slicer.load("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj");
			const char* Path = ("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj").c_str();
			remove(Path);
			slicer.save("postprocess\\" + mesh_target + "\\final_additive_blocks-" + to_string(cont_blocks) + ".obj");
			cont_blocks++;
		}
	}
	for (int i = 0; i < new_slicer_order.size() - 1; i++) {
		next_j = all_index_slice[new_slicer_order[i + 1].first].first;
		next_k = all_index_slice[new_slicer_order[i + 1].first].second;
		int j = all_index_slice[new_slicer_order[i].first].first;
		int k = all_index_slice[new_slicer_order[i].first].second;
		int jj = new_slicer_order[i].first;
		int kk = new_slicer_order[i].second;
		//for (int m = 0; m < all_slice_points[jj][kk].size(); m++) {
		//	out << "v " << all_slice_points[jj][kk][m].x << " " << all_slice_points[jj][kk][m].y << " " << all_slice_points[jj][kk][m].z << endl;
		//}
		cout << jj << " " << kk << endl;
		if (kk == all_slice_points[jj].size() - 1) {
			Slicer_2 slicer;
			slicer.load("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj");
			slicer.save("postprocess\\" + mesh_target + "\\final_additive_blocks-" + to_string(cont_blocks) + ".obj");
			const char* Path = ("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj").c_str();
			remove(Path);
			cont_blocks++;
			continue;
		}
		if (next_j != j || next_k != k) {
			vector<Vertex> temp_slice_points = all_slice_points[jj][kk];
			Slicer_2 slicer;
			Slicer_2 removed_triangles_from_slicer;
			Slicer_2 remain_triangles_from_slicer;
			slicer.load("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj");
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(all_slice_normal[jj][0].x, all_slice_normal[jj][0].y, all_slice_normal[jj][0].z);
			Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int m = 0; m < slicer.positions.size(); m++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = slicer.positions[m][0];
				temp_V(1, 0) = slicer.positions[m][1];
				temp_V(2, 0) = slicer.positions[m][2];
				temp_V = rotMatrix.inverse() * temp_V;
				slicer.positions[m][0] = temp_V(0, 0);
				slicer.positions[m][1] = temp_V(1, 0);
				slicer.positions[m][2] = temp_V(2, 0);
			}
			for (int m = 0; m < temp_slice_points.size(); m++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = temp_slice_points[m].x;
				temp_V(1, 0) = temp_slice_points[m].y;
				temp_V(2, 0) = temp_slice_points[m].z;
				temp_V = rotMatrix.inverse() * temp_V;
				temp_slice_points[m].x = temp_V(0, 0);
				temp_slice_points[m].y = temp_V(1, 0);
				temp_slice_points[m].z = temp_V(2, 0);
			}
			slicer.origin = { temp_slice_points[0].x,temp_slice_points[0].y ,temp_slice_points[0].z };
			slicer.normal = { 0,0,1 };
			slicer.cut();
			removed_triangles_from_slicer.positions = slicer.positions;
			remain_triangles_from_slicer.positions = slicer.positions;
			for (int m = 0; m < slicer.triangles.size(); m++) {
				for (int n = 0; n < 3; n++) {
					if (slicer.positions[slicer.triangles[m][n]][2] - slicer.origin[2] > 0.0000001) {
						removed_triangles_from_slicer.triangles.push_back(slicer.triangles[m]);
						slicer.triangles.erase(slicer.triangles.begin() + m);
						m--;
						break;
					}
				}
			}
			remain_triangles_from_slicer.triangles = slicer.triangles;
			remain_triangles_from_slicer.origin = slicer.origin;
			remain_triangles_from_slicer.normal = slicer.normal;
			Add_triangles_of_cutting_planes(remain_triangles_from_slicer);

			//旋转回来
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vectorBefore).toRotationMatrix();
			for (int m = 0; m < slicer.positions.size(); m++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = remain_triangles_from_slicer.positions[m][0];
				temp_V(1, 0) = remain_triangles_from_slicer.positions[m][1];
				temp_V(2, 0) = remain_triangles_from_slicer.positions[m][2];
				temp_V = rotMatrix.inverse() * temp_V;
				remain_triangles_from_slicer.positions[m][0] = temp_V(0, 0);
				remain_triangles_from_slicer.positions[m][1] = temp_V(1, 0);
				remain_triangles_from_slicer.positions[m][2] = temp_V(2, 0);
			}
			for (int m = 0; m < slicer.positions.size(); m++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = removed_triangles_from_slicer.positions[m][0];
				temp_V(1, 0) = removed_triangles_from_slicer.positions[m][1];
				temp_V(2, 0) = removed_triangles_from_slicer.positions[m][2];
				temp_V = rotMatrix.inverse() * temp_V;
				removed_triangles_from_slicer.positions[m][0] = temp_V(0, 0);
				removed_triangles_from_slicer.positions[m][1] = temp_V(1, 0);
				removed_triangles_from_slicer.positions[m][2] = temp_V(2, 0);
			}
			remain_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\final_additive_blocks-" + to_string(cont_blocks) + ".obj");
			removed_triangles_from_slicer.save("postprocess\\" + mesh_target + "\\additive_block-" + to_string(j) + "-" + to_string(k) + ".obj");
			cont_blocks++;
		}
	}

	cout << "Over!!" << endl;


}

bool PPS::Add_triangles_of_cutting_planes(Slicer_2& slicer)
{
	vector<vector<int>> cutting_plane_polygons;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(slicer.normal[0], slicer.normal[1], slicer.normal[2]);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	Slicer_2 temp_slicer = slicer;
	for (int j = 0; j < slicer.positions.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = slicer.positions[j][0];
		temp_V(1, 0) = slicer.positions[j][1];
		temp_V(2, 0) = slicer.positions[j][2];
		temp_V = rotMatrix.inverse() * temp_V;
		temp_slicer.positions[j][0] = temp_V(0, 0);
		temp_slicer.positions[j][1] = temp_V(1, 0);
		temp_slicer.positions[j][2] = temp_V(2, 0);
	}
	Eigen::MatrixXd temp_cutting_point;
	temp_cutting_point.resize(3, 1);
	temp_cutting_point(0, 0) = slicer.origin[0];
	temp_cutting_point(1, 0) = slicer.origin[1];
	temp_cutting_point(2, 0) = slicer.origin[2];
	temp_cutting_point = rotMatrix.inverse() * temp_cutting_point;
	vector<pair<int, int>> cutting_plane_edges;
	for (int j = 0; j < slicer.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if (abs(temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0)) < 0.0000001
				&& abs(temp_slicer.positions[slicer.triangles[j][(k + 1) % 3]][2] - temp_cutting_point(2, 0)) < 0.0000001)
				cutting_plane_edges.push_back(make_pair(slicer.triangles[j][k], slicer.triangles[j][(k + 1) % 3]));
		}
	}
	//Delete duplicate pairs in cutting_plane_edges
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		for (int k = j + 1; k < cutting_plane_edges.size(); k++) {
			if ((cutting_plane_edges[j].first == cutting_plane_edges[k].second && cutting_plane_edges[j].second == cutting_plane_edges[k].first)
				|| (cutting_plane_edges[j].first == cutting_plane_edges[k].first && cutting_plane_edges[j].second == cutting_plane_edges[k].second)) {
				cutting_plane_edges.erase(cutting_plane_edges.begin() + k);
				k--;
			}
		}
	}

	//get a polygon from the cutting_plane_edges,this polygon is made up of pairs in cutting_plane_edges,the adjacent points in a polygon are always in some pair
	vector<bool> flag_edge_been_search(cutting_plane_edges.size());
	for (int j = 0; j < cutting_plane_edges.size(); j++)
		flag_edge_been_search[j] = false;
	Point_3 temp_point(slicer.origin[0], slicer.origin[1], slicer.origin[2]);
	double min_distance = 100000000;
	int index_closest_point = -1;
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		double distance = std::sqrt(pow(temp_point.x() - slicer.positions[cutting_plane_edges[j].first][0], 2) + pow(temp_point.y() - slicer.positions[cutting_plane_edges[j].first][1], 2) + pow(temp_point.z() - slicer.positions[cutting_plane_edges[j].first][2], 2));
		if (min_distance > distance) {
			min_distance = distance;
			index_closest_point = j;
		}
	}
	vector<int> temp_vec;
	cutting_plane_polygons.push_back(temp_vec);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].first);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].second);
	flag_edge_been_search[index_closest_point] = true;
	int index_next_point = cutting_plane_edges[index_closest_point].second;
	while (true) {
		bool jud_find_next_point = false;
		for (int k = 0; k < cutting_plane_edges.size(); k++) {
			if (flag_edge_been_search[k] == true)
				continue;
			if (cutting_plane_edges[k].first == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			else if (cutting_plane_edges[k].second == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].first;
				jud_find_next_point = true;
				break;
			}
		}
		if (jud_find_next_point == false) {
			vector<int> temp_vec;
			cutting_plane_polygons.push_back(temp_vec);
			for (int k = 0; k < cutting_plane_edges.size(); k++) {
				if (flag_edge_been_search[k] == true)
					continue;
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			if (jud_find_next_point == false) {
				if (cutting_plane_polygons[cutting_plane_polygons.size() - 1].size() == 0)
					cutting_plane_polygons.erase(cutting_plane_polygons.end() - 1);
				break;
			}

		}
	}
	//return false if a cutting polygon is not closed
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
		if (cutting_plane_polygons[t][0] != cutting_plane_polygons[t][cutting_plane_polygons[t].size() - 1])
			return false;

	//add cutting planes to original mesh
	//Anticlockwise(cutting_plane_polygons, temp_slicer);
	//for (int i = 0; i < cutting_plane_polygons.size(); i++) {
	//	cutting_plane_polygons[i] = poufen(temp_slicer, cutting_plane_polygons[i], true);
	//	for (int j = 0; j < cutting_plane_polygons[i].size();) {
	//		TRiangle the_new_cutting_plane_triangle;
	//		the_new_cutting_plane_triangle[0] = cutting_plane_polygons[i][j]; j++;
	//		the_new_cutting_plane_triangle[1] = cutting_plane_polygons[i][j]; j++;
	//		the_new_cutting_plane_triangle[2] = cutting_plane_polygons[i][j]; j++;
	//		slicer.triangles.insert(slicer.triangles.end(), the_new_cutting_plane_triangle);  //添加接触面
	//	}
	//}

	//若不是洞上的边界，则删除
	for (int t = 0; t < cutting_plane_polygons.size(); t++) {
		int cont_num = 0;
		for (int j = 0; j < slicer.triangles.size(); j++)
			for (int k = 0; k < 3; k++) {
				if ((slicer.triangles[j][k] == cutting_plane_polygons[t][0] && slicer.triangles[j][(k + 1) % 3] == cutting_plane_polygons[t][1])
					|| (slicer.triangles[j][k] == cutting_plane_polygons[t][0] && slicer.triangles[j][(k + 2) % 3] == cutting_plane_polygons[t][1])) {
					cont_num++;
					break;
				}

			}
		if (cont_num > 1) {
			cutting_plane_polygons.erase(cutting_plane_polygons.begin() + t);
			t--;
		}
	}


	using Coord = double;
	using NN = uint32_t;
	using PPoint = std::array<Coord, 2>;
	std::vector<std::vector<PPoint>> polygon(1);
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
	{
		map<int, int> map_index_faces;
		map_index_faces.clear();
		polygon[0].clear();
		vector<Point_3> temp_vec;
		for (int j = 0; j < cutting_plane_polygons[t].size(); j++) {
			map_index_faces.insert({ j,cutting_plane_polygons[t][j] });
			PPoint temp_point;
			temp_point[0] = temp_slicer.positions[cutting_plane_polygons[t][j]][0];
			temp_point[1] = temp_slicer.positions[cutting_plane_polygons[t][j]][1];

			//防止由于共线导致的问题
			/*srand((unsigned)time(NULL));
			double temp_1 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[0] += temp_1;
			srand((unsigned)time(NULL));
			double temp_2 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[1] += temp_2;*/

			polygon[0].push_back(temp_point);
		}
		std::vector<NN> indices = mapbox::earcut<NN>(polygon);
		for (int j = 0; j < indices.size();) {
			TRiangle the_new_cutting_plane_triangle;
			the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
			slicer.triangles.insert(slicer.triangles.end(), the_new_cutting_plane_triangle);
		}
	}
	return true;
}
bool PPS::Add_triangles_of_cutting_planes_2(Slicer_2& slicer)
{
	vector<vector<int>> cutting_plane_polygons;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(slicer.normal[0], slicer.normal[1], slicer.normal[2]);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	Slicer_2 temp_slicer = slicer;
	for (int j = 0; j < slicer.positions.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = slicer.positions[j][0];
		temp_V(1, 0) = slicer.positions[j][1];
		temp_V(2, 0) = slicer.positions[j][2];
		temp_V = rotMatrix.inverse() * temp_V;
		temp_slicer.positions[j][0] = temp_V(0, 0);
		temp_slicer.positions[j][1] = temp_V(1, 0);
		temp_slicer.positions[j][2] = temp_V(2, 0);
	}
	Eigen::MatrixXd temp_cutting_point;
	temp_cutting_point.resize(3, 1);
	temp_cutting_point(0, 0) = slicer.origin[0];
	temp_cutting_point(1, 0) = slicer.origin[1];
	temp_cutting_point(2, 0) = slicer.origin[2];
	temp_cutting_point = rotMatrix.inverse() * temp_cutting_point;
	vector<pair<int, int>> cutting_plane_edges;
	for (int j = 0; j < slicer.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if (abs(temp_slicer.positions[slicer.triangles[j][k]][2] - temp_cutting_point(2, 0)) < 0.001
				&& abs(temp_slicer.positions[slicer.triangles[j][(k + 1) % 3]][2] - temp_cutting_point(2, 0)) < 0.001)
				cutting_plane_edges.push_back(make_pair(slicer.triangles[j][k], slicer.triangles[j][(k + 1) % 3]));
		}
	}
	//Delete duplicate pairs in cutting_plane_edges
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		for (int k = j + 1; k < cutting_plane_edges.size(); k++) {
			if ((cutting_plane_edges[j].first == cutting_plane_edges[k].second && cutting_plane_edges[j].second == cutting_plane_edges[k].first)
				|| (cutting_plane_edges[j].first == cutting_plane_edges[k].first && cutting_plane_edges[j].second == cutting_plane_edges[k].second)) {
				cutting_plane_edges.erase(cutting_plane_edges.begin() + k);
				k--;
			}
		}
	}

	//get a polygon from the cutting_plane_edges,this polygon is made up of pairs in cutting_plane_edges,the adjacent points in a polygon are always in some pair
	vector<bool> flag_edge_been_search(cutting_plane_edges.size());
	for (int j = 0; j < cutting_plane_edges.size(); j++)
		flag_edge_been_search[j] = false;
	Point_3 temp_point(slicer.origin[0], slicer.origin[1], slicer.origin[2]);
	double min_distance = 100000000;
	int index_closest_point = -1;
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		double distance = std::sqrt(pow(temp_point.x() - slicer.positions[cutting_plane_edges[j].first][0], 2) + pow(temp_point.y() - slicer.positions[cutting_plane_edges[j].first][1], 2) + pow(temp_point.z() - slicer.positions[cutting_plane_edges[j].first][2], 2));
		if (min_distance > distance) {
			min_distance = distance;
			index_closest_point = j;
		}
	}
	vector<int> temp_vec;
	cutting_plane_polygons.push_back(temp_vec);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].first);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].second);
	flag_edge_been_search[index_closest_point] = true;
	int index_next_point = cutting_plane_edges[index_closest_point].second;
	while (true) {
		bool jud_find_next_point = false;
		for (int k = 0; k < cutting_plane_edges.size(); k++) {
			if (flag_edge_been_search[k] == true)
				continue;
			if (cutting_plane_edges[k].first == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			else if (cutting_plane_edges[k].second == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].first;
				jud_find_next_point = true;
				break;
			}
		}
		if (jud_find_next_point == false) {
			vector<int> temp_vec;
			cutting_plane_polygons.push_back(temp_vec);
			for (int k = 0; k < cutting_plane_edges.size(); k++) {
				if (flag_edge_been_search[k] == true)
					continue;
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			if (jud_find_next_point == false) {
				if (cutting_plane_polygons[cutting_plane_polygons.size() - 1].size() == 0)
					cutting_plane_polygons.erase(cutting_plane_polygons.end() - 1);
				break;
			}

		}
	}
	//return false if a cutting polygon is not closed
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
		if (cutting_plane_polygons[t][0] != cutting_plane_polygons[t][cutting_plane_polygons[t].size() - 1])
			return false;

	//add cutting planes to original mesh
	//Anticlockwise(cutting_plane_polygons, temp_slicer);
	//for (int i = 0; i < cutting_plane_polygons.size(); i++) {
	//	cutting_plane_polygons[i] = poufen(temp_slicer, cutting_plane_polygons[i], true);
	//	for (int j = 0; j < cutting_plane_polygons[i].size();) {
	//		TRiangle the_new_cutting_plane_triangle;
	//		the_new_cutting_plane_triangle[0] = cutting_plane_polygons[i][j]; j++;
	//		the_new_cutting_plane_triangle[1] = cutting_plane_polygons[i][j]; j++;
	//		the_new_cutting_plane_triangle[2] = cutting_plane_polygons[i][j]; j++;
	//		slicer.triangles.insert(slicer.triangles.end(), the_new_cutting_plane_triangle);  //添加接触面
	//	}
	//}

	//若不是洞上的边界，则删除
	for (int t = 0; t < cutting_plane_polygons.size(); t++) {
		int cont_num = 0;
		for (int j = 0; j < slicer.triangles.size(); j++)
			for (int k = 0; k < 3; k++) {
				if ((slicer.triangles[j][k] == cutting_plane_polygons[t][0] && slicer.triangles[j][(k + 1) % 3] == cutting_plane_polygons[t][1])
					|| (slicer.triangles[j][k] == cutting_plane_polygons[t][0] && slicer.triangles[j][(k + 2) % 3] == cutting_plane_polygons[t][1])) {
					cont_num++;
					break;
				}

			}
		if (cont_num > 1) {
			cutting_plane_polygons.erase(cutting_plane_polygons.begin() + t);
			t--;
		}
	}


	using Coord = double;
	using NN = uint32_t;
	using PPoint = std::array<Coord, 2>;
	std::vector<std::vector<PPoint>> polygon(1);
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
	{
		map<int, int> map_index_faces;
		map_index_faces.clear();
		polygon[0].clear();
		vector<Point_3> temp_vec;
		for (int j = 0; j < cutting_plane_polygons[t].size(); j++) {
			map_index_faces.insert({ j,cutting_plane_polygons[t][j] });
			PPoint temp_point;
			temp_point[0] = temp_slicer.positions[cutting_plane_polygons[t][j]][0];
			temp_point[1] = temp_slicer.positions[cutting_plane_polygons[t][j]][1];

			//防止由于共线导致的问题
			/*srand((unsigned)time(NULL));
			double temp_1 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[0] += temp_1;
			srand((unsigned)time(NULL));
			double temp_2 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[1] += temp_2;*/

			polygon[0].push_back(temp_point);
		}
		std::vector<NN> indices = mapbox::earcut<NN>(polygon);
		for (int j = 0; j < indices.size();) {
			TRiangle the_new_cutting_plane_triangle;
			the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
			slicer.triangles.insert(slicer.triangles.end(), the_new_cutting_plane_triangle);
		}
	}
	return true;
}

bool PPS::Check_cutting_bottom(Slicer_2& slicer)
{
	vector<vector<int>> cutting_plane_polygons;
	vector<pair<int, int>> cutting_plane_edges;
	for (int j = 0; j < slicer.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if (abs(slicer.positions[slicer.triangles[j][k]][2] - slicer.origin[2]) < 0.00001
				&& abs(slicer.positions[slicer.triangles[j][(k + 1) % 3]][2] - slicer.origin[2]) < 0.00001)
				cutting_plane_edges.push_back(make_pair(slicer.triangles[j][k], slicer.triangles[j][(k + 1) % 3]));
		}
	}
	//Delete duplicate pairs in cutting_plane_edges
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		for (int k = j + 1; k < cutting_plane_edges.size(); k++) {
			if ((cutting_plane_edges[j].first == cutting_plane_edges[k].second && cutting_plane_edges[j].second == cutting_plane_edges[k].first)
				|| (cutting_plane_edges[j].first == cutting_plane_edges[k].first && cutting_plane_edges[j].second == cutting_plane_edges[k].second)) {
				cutting_plane_edges.erase(cutting_plane_edges.begin() + k);
				k--;
			}
		}
	}

	//get a polygon from the cutting_plane_edges,this polygon is made up of pairs in cutting_plane_edges,the adjacent points in a polygon are always in some pair
	vector<bool> flag_edge_been_search(cutting_plane_edges.size());
	for (int j = 0; j < cutting_plane_edges.size(); j++)
		flag_edge_been_search[j] = false;
	Point_3 temp_point(slicer.origin[0], slicer.origin[1], slicer.origin[2]);
	double min_distance = 100000000;
	int index_closest_point = -1;
	for (int j = 0; j < cutting_plane_edges.size(); j++) {
		double distance = std::sqrt(pow(temp_point.x() - slicer.positions[cutting_plane_edges[j].first][0], 2) + pow(temp_point.y() - slicer.positions[cutting_plane_edges[j].first][1], 2) + pow(temp_point.z() - slicer.positions[cutting_plane_edges[j].first][2], 2));
		if (min_distance > distance) {
			min_distance = distance;
			index_closest_point = j;
		}
	}
	vector<int> temp_vec;
	cutting_plane_polygons.push_back(temp_vec);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].first);
	cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[index_closest_point].second);
	flag_edge_been_search[index_closest_point] = true;
	int index_next_point = cutting_plane_edges[index_closest_point].second;
	while (true) {
		bool jud_find_next_point = false;
		for (int k = 0; k < cutting_plane_edges.size(); k++) {
			if (flag_edge_been_search[k] == true)
				continue;
			if (cutting_plane_edges[k].first == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			else if (cutting_plane_edges[k].second == index_next_point) {
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].first;
				jud_find_next_point = true;
				break;
			}
		}
		if (jud_find_next_point == false) {
			vector<int> temp_vec;
			cutting_plane_polygons.push_back(temp_vec);
			for (int k = 0; k < cutting_plane_edges.size(); k++) {
				if (flag_edge_been_search[k] == true)
					continue;
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].first);
				cutting_plane_polygons[cutting_plane_polygons.size() - 1].push_back(cutting_plane_edges[k].second);
				flag_edge_been_search[k] = true;
				index_next_point = cutting_plane_edges[k].second;
				jud_find_next_point = true;
				break;
			}
			if (jud_find_next_point == false) {
				if (cutting_plane_polygons[cutting_plane_polygons.size() - 1].size() == 0)
					cutting_plane_polygons.erase(cutting_plane_polygons.end() - 1);
				break;
			}

		}
	}
	//return false if a cutting polygon is not closed
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
		if (cutting_plane_polygons[t][0] != cutting_plane_polygons[t][cutting_plane_polygons[t].size() - 1])
			return false;

	using Coord = double;
	using NN = uint32_t;
	using PPoint = std::array<Coord, 2>;
	std::vector<std::vector<PPoint>> polygon(1);
	for (int t = 0; t < cutting_plane_polygons.size(); t++)
	{
		map<int, int> map_index_faces;
		map_index_faces.clear();
		polygon[0].clear();
		vector<Point_3> temp_vec;
		for (int j = 0; j < cutting_plane_polygons[t].size(); j++) {
			map_index_faces.insert({ j,cutting_plane_polygons[t][j] });
			PPoint temp_point;
			temp_point[0] = slicer.positions[cutting_plane_polygons[t][j]][0];
			temp_point[1] = slicer.positions[cutting_plane_polygons[t][j]][1];

			//防止由于共线导致的问题
			/*srand((unsigned)time(NULL));
			double temp_1 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[0] += temp_1;
			srand((unsigned)time(NULL));
			double temp_2 = static_cast<double>(rand()) / RAND_MAX * 0.00000001;
			temp_point[1] += temp_2;*/

			polygon[0].push_back(temp_point);
		}
		std::vector<NN> indices = mapbox::earcut<NN>(polygon);
		for (int j = 0; j < indices.size();) {
			TRiangle the_new_cutting_plane_triangle;
			the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
			the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
			slicer.triangles.insert(slicer.triangles.end(), the_new_cutting_plane_triangle);
		}
	}
	return true;
}

void PPS::OrientationSamplePoints()
{
	sample_points_sphere.clear();
	int number_points = 50;  //50  
	sample_points_sphere.push_back(Point_3(0, 0, 1));
	double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
	double offset = 2.0 / number_points;
	for (int i = 0; i < number_points; ++i) {
		double y = i * offset - 1.0 + (offset / 2.0);
		double radiusXY = sqrt(1.0 - y * y);
		double theta = goldenAngle * i;
		double x = cos(theta) * radiusXY;
		double z = sin(theta) * radiusXY;
		if (z > 0)
			sample_points_sphere.push_back({ x * 1, y * 1, z * 1 });
	}
}

double PPS::Calculate_area_need_support(Slicer_2& slicer, double& sum_area)
{
	double hang_in_area = 0;
	vector<Vector3> vertices_in_face;
	for (int i = 0; i < slicer.triangles.size(); i++) {
		vertices_in_face.clear();
		for (int j = 0; j < 3; j++) {
			vertices_in_face.push_back(Vector3(slicer.positions[slicer.triangles[i][j]][0], slicer.positions[slicer.triangles[i][j]][1], slicer.positions[slicer.triangles[i][j]][2]));
		}
		Vector3 v1 = vertices_in_face[0];
		Vector3 v2 = vertices_in_face[1];
		Vector3 v3 = vertices_in_face[2];

		double na = (v2.m_y - v1.m_y) * (v3.m_z - v1.m_z) - (v2.m_z - v1.m_z) * (v3.m_y - v1.m_y);
		double nb = (v2.m_z - v1.m_z) * (v3.m_x - v1.m_x) - (v2.m_x - v1.m_x) * (v3.m_z - v1.m_z);
		double nc = (v2.m_x - v1.m_x) * (v3.m_y - v1.m_y) - (v2.m_y - v1.m_y) * (v3.m_x - v1.m_x);
		Vector3 face_normal(na, nb, nc);
		face_normal.Normalized();
		Vector3 base_normal(0, 0, 1);
		base_normal.Normalized();
		bool jud_self_support;
		if (Dot(face_normal, base_normal) + sin(PI / 3.6) >= 0)
			jud_self_support = true;
		else
			jud_self_support = false;
		Eigen::Vector3d vec1, vec2;
		vec1.x() = v2.m_x - v1.m_x;
		vec1.y() = v2.m_y - v1.m_y;
		vec1.z() = v2.m_z - v1.m_z;
		vec2.x() = v3.m_x - v1.m_x;
		vec2.y() = v3.m_y - v1.m_y;
		vec2.z() = v3.m_z - v1.m_z;
		if (jud_self_support == false)
		{
			hang_in_area += vec1.cross(vec2).norm() / 2;
		}
		sum_area += vec1.cross(vec2).norm() / 2;
	}

	return hang_in_area;
}

double PPS::Calculate_area_need_support_2(Slicer_2& slicer, double& sum_area, bool& jud_flat_area)
{
	double hang_in_area = 0;
	vector<Vector3> vertices_in_face;
	for (int i = 0; i < slicer.triangles.size(); i++) {
		vertices_in_face.clear();
		for (int j = 0; j < 3; j++) {
			vertices_in_face.push_back(Vector3(slicer.positions[slicer.triangles[i][j]][0], slicer.positions[slicer.triangles[i][j]][1], slicer.positions[slicer.triangles[i][j]][2]));
		}
		Vector3 v1 = vertices_in_face[0];
		Vector3 v2 = vertices_in_face[1];
		Vector3 v3 = vertices_in_face[2];

		double na = (v2.m_y - v1.m_y) * (v3.m_z - v1.m_z) - (v2.m_z - v1.m_z) * (v3.m_y - v1.m_y);
		double nb = (v2.m_z - v1.m_z) * (v3.m_x - v1.m_x) - (v2.m_x - v1.m_x) * (v3.m_z - v1.m_z);
		double nc = (v2.m_x - v1.m_x) * (v3.m_y - v1.m_y) - (v2.m_y - v1.m_y) * (v3.m_x - v1.m_x);
		Vector3 face_normal(na, nb, nc);
		face_normal.Normalized();
		Vector3 base_normal(0, 0, 1);
		base_normal.Normalized();
		bool jud_self_support;
		if (Dot(face_normal, base_normal) + sin(PI / 2.8) < 0) {    //2.8
			jud_flat_area = true;
			break;
		}
		if (Dot(face_normal, base_normal) + sin(PI / 3.6) >= 0)
			jud_self_support = true;
		else
			jud_self_support = false;
		Eigen::Vector3d vec1, vec2;
		vec1.x() = v2.m_x - v1.m_x;
		vec1.y() = v2.m_y - v1.m_y;
		vec1.z() = v2.m_z - v1.m_z;
		vec2.x() = v3.m_x - v1.m_x;
		vec2.y() = v3.m_y - v1.m_y;
		vec2.z() = v3.m_z - v1.m_z;
		if (jud_self_support == false)
		{
			hang_in_area += vec1.cross(vec2).norm() / 2;
		}
		sum_area += vec1.cross(vec2).norm() / 2;
	}

	return hang_in_area;
}