#include"ReFab.h"

//boolean:difference
void ReFab::Intersection_two_mesh(std::string& input1, std::string& input2, std::string output,bool type)
{
	if (type == 0) {
		//读入mesh1和mesh2
		Mesh mesh1, mesh2, inter;
		if (!PMP::IO::read_polygon_mesh(rootPath + input1, mesh1) || !PMP::IO::read_polygon_mesh(rootPath + input2, mesh2))
		{
			std::cerr << "Invalid input." << std::endl;
			return;
		}
		Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
		Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
		Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
		//对于两个mesh求交，相交的mesh保存到inter
		if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
			params::vertex_point_map(mesh1_vpm),
			params::vertex_point_map(mesh2_vpm),
			params::vertex_point_map(inter_vpm)))
		{
			std::cout << "Intersection were successfully computed\n";
			CGAL::IO::write_polygon_mesh(rootPath + output, inter, CGAL::parameters::stream_precision(17));
		}
		else {
			std::cout << "Intersection could not be computed\n";
			return;
		}
	}
	

	else {
		Eigen::MatrixXd VA, VB, VC;
		Eigen::VectorXi J, I;
		Eigen::MatrixXi FA, FB, FC;
		igl::MeshBooleanType boolean_type(
			igl::MESH_BOOLEAN_TYPE_INTERSECT);
		igl::readOFF(rootPath + input1, VA, FA);
		igl::readOFF(rootPath + input2, VB, FB);
		igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
		igl::writeOFF(rootPath + output, VC, FC);
	}
	
}

//boolean:intersection
void ReFab::Difference_two_mesh(std::string& input1, std::string& input2, std::string output, int type)
{
	//读入mesh1和mesh2
	if (type == 0) {
		Mesh mesh1, mesh2, subtract;
		if (!PMP::IO::read_polygon_mesh(rootPath + input1, mesh1) || !PMP::IO::read_polygon_mesh(rootPath + input2, mesh2))
		{
			std::cerr << "Invalid input." << std::endl;
			return;
		}
		//对两个mesh做差
		//create a property on edges to indicate whether they are constrained
		//Mesh::Property_map<edge_descriptor, bool> is_constrained_map = subtract.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
		Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		Exact_point_map inter_exact_points = subtract.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
		
		

		Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
		Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
		Exact_vertex_point_map inter_vpm(inter_exact_points, subtract);

		// update mesh1 to contain the mesh bounding the difference
		// of the two input volumes.
		if (PMP::corefine_and_compute_difference(mesh1, mesh2, subtract,
			params::vertex_point_map(mesh1_vpm),
			params::vertex_point_map(mesh2_vpm),
			params::vertex_point_map(inter_vpm)))
		{
			std::cout << "Difference1 was successfully computed\n";
			CGAL::IO::write_polygon_mesh(rootPath + output, subtract, CGAL::parameters::stream_precision(17));
		}
		else
		{
			std::cout << "Difference1 could not be computed\n";
			return;
		}
	}
	
	else {
		Eigen::MatrixXd VA, VB, VC;
		Eigen::VectorXi J, I;
		Eigen::MatrixXi FA, FB, FC;
		igl::MeshBooleanType boolean_type(
			igl::MESH_BOOLEAN_TYPE_MINUS);
		igl::readOFF(rootPath + input1, VA, FA);
		igl::readOFF(rootPath + input2, VB, FB);
		igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
		igl::writeOFF(rootPath + output, VC, FC);
	}
}

void ReFab::DetermineBase(string path)
{
	//********用的是初始形状的基底*********//

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	using namespace Eigen;
	igl::read_triangle_mesh(path, V, F);
	//calculate the minimum z value of all the triangles
	double min_z = 100000;
	for (int i = 0; i < F.rows(); i++) {
		int index0 = F(i, 0);
		int index1 = F(i, 1);
		int index2 = F(i, 2);
		double z0 = V.row(index0).z();
		double z1 = V.row(index1).z();
		double z2 = V.row(index2).z();
		if (z0 < min_z) min_z = z0;
		if (z1 < min_z) min_z = z1;
		if (z2 < min_z) min_z = z2;
	}

	//计算法向
	for (int i = 0; i < F.rows(); i++) {
		int index0 = F(i, 0);
		int index1 = F(i, 1);
		int index2 = F(i, 2);
		double na = (V.row(index1).y() - V.row(index0).y()) * (V.row(index2).z() - V.row(index0).z()) - (V.row(index1).z() - V.row(index0).z()) * (V.row(index2).y() - V.row(index0).y());
		double nb = (V.row(index1).z() - V.row(index0).z()) * (V.row(index2).x() - V.row(index0).x()) - (V.row(index1).x() - V.row(index0).x()) * (V.row(index2).z() - V.row(index0).z());
		double nc = (V.row(index1).x() - V.row(index0).x()) * (V.row(index2).y() - V.row(index0).y()) - (V.row(index1).y() - V.row(index0).y()) * (V.row(index2).x() - V.row(index0).x());
		Vector3 vn(na, nb, nc);
		vn.Normalized();
		if(abs(V.row(F(i, 0)).z() - min_z) < 0.001 && abs(V.row(F(i, 1)).z() - min_z) < 0.001 && abs(V.row(F(i, 2)).z() - min_z) < 0.001)
			if (abs(vn.m_x)< 0.001  && abs(vn.m_y) < 0.001 && abs(vn.m_z + 1) < 0.001) {
				vertices_base.push_back(Vec3(V.row(index0).x(), V.row(index0).y(), V.row(index0).z()));
				vertices_base.push_back(Vec3(V.row(index1).x(), V.row(index1).y(), V.row(index1).z()));
				vertices_base.push_back(Vec3(V.row(index2).x(), V.row(index2).y(), V.row(index2).z()));
				index_base_of_Di.push_back(i);
			}
	}
}

//boolean:union
void Union_two_mesh(std::string& input1, std::string& input2, std::string output)
{
	//ToDO
}

void ReFab::CalculateHangInArea(string input_DI)
{
	Mesh mesh,mesh2;
	PMP::IO::read_polygon_mesh(input_DI, mesh);
	Eigen::MatrixXd V;
	Eigen::MatrixXi F, F_new, DI_F_new;
	igl::readOFF(input_DI, V, F);

	double min_z = 100000;
	face_descriptor min_face;
	vector<bool> flag_face_been_search(mesh.number_of_faces());
	for (int i = 0; i < mesh.number_of_faces(); i++)
		flag_face_been_search[i] = false;
	typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
	const double bound = std::cos(0.75 * CGAL_PI);
	for (face_descriptor fd : faces(mesh))
	{
		for (Mesh::Vertex_index vd : mesh.vertices_around_face(mesh.halfedge(fd)))
		{
			Point_3 p = mesh.point(vd);
			if (p.z() < min_z) {
				min_z = p.z();
				min_face = fd;
			}
		}
	}
	std::vector<face_descriptor> cc;
	cc.clear();
	PMP::connected_component(min_face,
		mesh,
		std::back_inserter(cc));


	bool* flag_hang_in_face = new bool[mesh.num_faces()];
	for (int i = 0; i < mesh.num_faces(); i++)
		flag_hang_in_face[i] = true;
	DI_F_new.resize(cc.size(), 3);
	for (int i = 0; i < cc.size(); i++) {
		flag_hang_in_face[cc[i].idx()] = false;
		DI_F_new.row(i) = F.row(cc[i].idx());
	}
	F_new.resize(mesh.num_faces() - cc.size(), 3);
	int index_F_new = 0;
	for (int i = 0; i < mesh.num_faces(); i++) {
		if (flag_hang_in_face[i] == true) {
			F_new.row(index_F_new) = F.row(i);
			index_F_new++;
		}
	}

	igl::writeOFF(rootPath + "output\\"+mesh_target+"\\hang_in_area.off", V, F_new);
	igl::writeOFF(rootPath + "output\\"+mesh_target+"\\D_I.off", V, DI_F_new);
	//PMP::IO::read_polygon_mesh(rootPath + "output\\D_I.off", mesh2);
	//CGAL::IO::write_polygon_mesh(rootPath + "output\\D_I.off", mesh2, CGAL::parameters::stream_precision(17));
}

void ReFab::Modify_DS_and_DA(string input_DS, string input_DA, string input_hang_in_area)
{
	Mesh mesh1, mesh2, mesh3, mesh4, union_mesh, union_mesh_2;
	if (!PMP::IO::read_polygon_mesh(rootPath + input_DS, mesh1) || !PMP::IO::read_polygon_mesh(rootPath + input_hang_in_area, mesh2))
	{
		std::cerr << "Invalid input." << std::endl;
		return;
	}
	Mesh::Property_map<edge_descriptor, bool> is_constrained_map = union_mesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
	if (PMP::corefine_and_compute_union(mesh1, mesh2, union_mesh,
		params::default_values(), // default parameters for mesh1
		params::default_values(), // default parameters for mesh2
		params::edge_is_constrained_map(is_constrained_map)))
	{
		CGAL::IO::write_polygon_mesh(rootPath + input_DS, union_mesh, CGAL::parameters::stream_precision(17));
	}
	else
	{
		std::cout << "Difference1 could not be computed\n";
		return;
	}

	if (!PMP::IO::read_polygon_mesh(rootPath + input_DA, mesh3) || !PMP::IO::read_polygon_mesh(rootPath + input_hang_in_area, mesh4))
	{
		std::cerr << "Invalid input." << std::endl;
		return;
	}
	Mesh::Property_map<edge_descriptor, bool> is_constrained_map_2 = union_mesh_2.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
	if (PMP::corefine_and_compute_union(mesh3, mesh4, union_mesh_2,
		params::default_values(), // default parameters for mesh1
		params::default_values(), // default parameters for mesh2
		params::edge_is_constrained_map(is_constrained_map_2)))
	{
		CGAL::IO::write_polygon_mesh(rootPath + input_DA, union_mesh_2, CGAL::parameters::stream_precision(17));
	}
	else
	{
		std::cout << "Difference1 could not be computed\n";
		return;
	}
}

//对模型进行求交，并进行分割，返回分割后的模型
void ReFab::Intersection_and_Split(std::string& input_s, std::string& input_a, std::string output_s, std::string output_a, std::string output_inter,bool type)
{

	//对mesh_a和mesh_s求交
	Intersection_two_mesh(input_s, input_a, output_inter,type);
	CalculateHangInArea(rootPath + output_inter);

	Mesh mesh_a, mesh_s, inter;
	if (!PMP::IO::read_polygon_mesh(rootPath + input_s, mesh_s) || !PMP::IO::read_polygon_mesh(rootPath + input_a, mesh_a) || !PMP::IO::read_polygon_mesh(rootPath + output_inter, inter))
	{
		std::cerr << "Invalid input." << std::endl;
		return;
	}
	volume_of_target_model = CGAL::Polygon_mesh_processing::volume(mesh_a);
	record_data[0] = (CGAL::Polygon_mesh_processing::volume(inter) / volume_of_target_model);

	std::cout << "Spliting\n";
	//判断，如果inter的面在mesh1内部，则将inter保存到output1中，将mesh1-inter保存到output2中
	//对inter中的每一个面都进行判断
	//接触面分为增材接触面和减材接触面
	Mesh interface_s, interface_a;

	//建立一个新的映射，对于interface_s中的每一个点，都有一个对应的点在inter中
	std::unordered_map<Mesh::Vertex_index, Mesh::Vertex_index> interface_s_map, interface_a_map;

	//创建一个内外判别器，用于判断点是否在mesh1内部
	CGAL::Side_of_triangle_mesh<Mesh, K> inside_s(mesh_s), inside_a(mesh_a);
	//遍历inter中的每一个面
	double min_z_of_DI = 100000;
	for (Mesh::Face_index face : faces(inter)) {
		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			if (inter.point(vd).z() < min_z_of_DI)
				min_z_of_DI = inter.point(vd).z();
		}
	}

	//遍历inter中的每一个面
	for (Mesh::Face_index face : faces(inter))
	{
		//就只计算一个面的终点是不是在mesh_a内部

		double x = 0.0, y = 0.0, z = 0.0;
		std::vector<Point_3> PointOfFace;

		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			Point_3 p = inter.point(vd);
			PointOfFace.push_back(p);
			//将xyz累加起来
			x += p.x();
			y += p.y();
			z += p.z();
		}
		//计算face的法向
		K::Vector_3 normal = CGAL::normal(PointOfFace[0], PointOfFace[1], PointOfFace[2]);
		K::Vector_3 center(x / 3.0, y / 3.0, z / 3.0);
		center = center + 0.01 * normal;


		Point_3 centerP(center.x(), center.y(), center.z());
		//判断中心点是否在mesh_a内部
		if (inside_a(centerP) == CGAL::ON_CONVEX_SIDE) {//说明这个面落在mesh_a内部
			//存储该面对应的点的索引
			std::vector<Mesh::Vertex_index> new_vertices;
			//将组成面的点保存到interface_a中
			for (Mesh::Vertex_index vertex_index : vertices_around_face(inter.halfedge(face), inter)) {
				//如果顶点还没有添加到interface_s，则添加它
				if (interface_a_map.find(vertex_index) == interface_a_map.end()) {
					Point_3 p = inter.point(vertex_index);
					Mesh::Vertex_index new_vertex = interface_a.add_vertex(p);
					interface_a_map[vertex_index] = new_vertex;
				}
				new_vertices.push_back(interface_a_map[vertex_index]);
			}

			//添加新的面到interface_s中,对于定点数不同的面分别处理
			if (new_vertices.size() == 3)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2]);
			else if (new_vertices.size() == 4)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2], new_vertices[3]);
			else
				std::cerr << "Error: Unable to add face to interface_a." << std::endl;

		}
		else {
			//存储该面对应的点的索引
			std::vector<Mesh::Vertex_index> new_vertices;
			//将组成面的点保存到interface_s中
			for (Mesh::Vertex_index vertex_index : vertices_around_face(inter.halfedge(face), inter)) {
				//如果顶点还没有添加到interface_s，则添加它
				if (interface_s_map.find(vertex_index) == interface_s_map.end()) {
					Point_3 p = inter.point(vertex_index);
					Mesh::Vertex_index new_vertex = interface_s.add_vertex(p);
					interface_s_map[vertex_index] = new_vertex;
				}
				new_vertices.push_back(interface_s_map[vertex_index]);
			}

			//添加新的面到interface_s中,对于定点数不同的面分别处理
			if (new_vertices.size() == 3)
				Mesh::Face_index new_face = interface_s.add_face(new_vertices[0], new_vertices[1], new_vertices[2]);
			else if (new_vertices.size() == 4)
				Mesh::Face_index new_face = interface_s.add_face(new_vertices[0], new_vertices[1], new_vertices[2], new_vertices[3]);
			else
				std::cerr << "Error: Unable to add face to interface_s." << std::endl;
		}
	}

	//将interface_s保存到output1中
	std::cout << "Split finished\n";
	CGAL::IO::write_polygon_mesh(rootPath + output_s, interface_s, CGAL::parameters::stream_precision(17));
	CGAL::IO::write_polygon_mesh(rootPath + output_a, interface_a, CGAL::parameters::stream_precision(17));
}


//找到给定三角形网格（Polyhedron 类型）中所有顶点的最大坐标值
double max_coordinate(const Mesh& mesh)
{
	double max_coord = -std::numeric_limits<double>::infinity();
	for (Mesh::Vertex_index v : vertices(mesh))
	{
		Point_3 p = mesh.point(v);
		max_coord = (std::max)(max_coord, CGAL::to_double(p.x()));
		max_coord = (std::max)(max_coord, CGAL::to_double(p.y()));
		max_coord = (std::max)(max_coord, CGAL::to_double(p.z()));
	}
	return max_coord;
}


// 辅助函数：判断面是否属于指定的网格
// 通过face的顶点坐标的对应关系来判断
bool is_face_in_mesh(const Mesh& mesh, const Mesh::Face_index& face)
{
	//创建一个inside，用于判断点是否在mesh内部
	CGAL::Side_of_triangle_mesh<Mesh, K> inside(mesh);

	std::vector<Point_3> PointOfFace;
	for (Mesh::Vertex_index vd : mesh.vertices_around_face(mesh.halfedge(face)))
		PointOfFace.push_back(mesh.point(vd));

	//判断面的三个顶点是否在mesh内部
	//如果三个点都在mesh内部，则该面在mesh内部，返回true
	//否则返回false
	for (auto p : PointOfFace)
	{
		CGAL::Bounded_side res = inside(p);
		//有一个点在mesh的边界上,或者在边界外，则该面不在mesh内部
		if (res == CGAL::ON_BOUNDARY || res == CGAL::ON_UNBOUNDED_SIDE) { return false; }
	}

	return true;
}


/********************************************************
* Function name: CGAL_3D_Mesh_Regular_Sampling_C1
* Description: This function performs regular grid-based sampling within a three-dimensional polyhedral mesh.
*				It generates a set of sampling points at regular intervals within the specified region.
* Parameters:
* @outside_path: A string representing the path or filename of the mesh to be sampled.
* @d: percentage value of the length of the diagonal of the bounding box.
* @sampling_points: An output parameter, a vector of three-dimensional points representing the generated sampling points within the specified region.
* Return: void
**********************************************************/
void ReFab::Mesh_Regular_Sampling_C1(std::string& outside_path, const double& d)
{
	sampling_points_in_D_I.clear();


	// 读入outside_mesh
	Mesh outside_mesh;
	if (!PMP::IO::read_polygon_mesh(rootPath + outside_path, outside_mesh))
	{
		std::cerr << "Invalid input." << std::endl;
		return;
	}
	// 计算边界框 
	CGAL::Bbox_3 bbox = PMP::bbox(outside_mesh);
	Point_3 out_minC(bbox.xmin(), bbox.ymin(), bbox.zmin()), out_maxC(bbox.xmax(), bbox.ymax(), bbox.zmax());

	//std::cout << "Bounding Box Min: (" << out_minC[0] << ", " << out_minC[1] << ", " << out_minC[2] << ")" << std::endl;
	//std::cout << "Bounding Box Max: (" << out_maxC[0] << ", " << out_maxC[1] << ", " << out_maxC[2] << ")" << std::endl;
	std::cout << "regular sampling start" << std::endl;
	double diagonal_length = sqrt(CGAL::squared_distance(out_minC, out_maxC));
	//直接用d做采样步长
	//minimal_d = diagonal_length * d;
	minimal_d = 2;  //4 //2

	CGAL::Side_of_triangle_mesh<Mesh, K> inside(outside_mesh);
	//沿着xyz轴方向，以minimal_d为步长，生成采样点
	double x(out_minC[0]);
	while (x < out_maxC[0])
	{
		double y(out_minC[1]);
		while (y < out_maxC[1])
		{
			double z(out_minC[2]);
			while (z < out_maxC[2])
			{
				Point_3 p(x, y, z);
				CGAL::Bounded_side res = inside(p);
				if (res == CGAL::ON_BOUNDED_SIDE)
				{
					sampling_points_in_D_I.push_back(p);
					if (sampling_points_in_D_I.size() % 10000 == 0) std::cout << sampling_points_in_D_I.size() << " ";
				}
				z += minimal_d;
			}
			y += minimal_d;
		}
		x += minimal_d;
	}

	///////////////////////////////////暂时用于测试插值,之后注释////////////////////////////
	//float minimal_d2 = 2;  //4
	////沿着xyz轴方向，以minimal_d为步长，生成采样点
	//x = out_minC[0];
	//while (x < out_maxC[0])
	//{
	//	double y(out_minC[1]);
	//	while (y < out_maxC[1])
	//	{
	//		double z(out_minC[2]);
	//		while (z < out_maxC[2])
	//		{
	//			Point_3 p(x, y, z);
	//			CGAL::Bounded_side res = inside(p);
	//			if (res == CGAL::ON_BOUNDED_SIDE)
	//			{
	//				Eigen::Vector3d temp(p.x(), p.y(), p.z());
	//				test_sampling_points.push_back(temp);
	//			}
	//			z += minimal_d2;
	//		}
	//		y += minimal_d2;
	//	}
	//	x += minimal_d2;
	//}
	//////////////////////////////////////////////////////////////////////////////

	ofstream ofile(rootPath + "output\\"+mesh_target+"\\D_I_sampling_points.obj");
	for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
		ofile << "v " << sampling_points_in_D_I[i].x() << " " << sampling_points_in_D_I[i].y() << " " << sampling_points_in_D_I[i].z() << std::endl;
	}
	std::cerr << std::endl;
}

void ReFab::Get_sampling_points_in_DS_and_DA(std::string DS_file, std::string DA_file)
{
	ifstream DS_points(rootPath+DS_file);
	double x, y, z,nx,ny,nz;
	while (!DS_points.eof()) {
		DS_points >> x >> y >> z;
		Point_3 temp_point = Point_3(x, y, z);
		sampling_points_in_D_S.push_back(temp_point);
		DS_points >> nx >> ny >> nz;
		Point_3 temp_normal = Point_3(nx, ny, nz);
		normal_D_S.push_back(temp_normal);
	}
	ifstream DA_points(rootPath+DA_file);
	while (!DA_points.eof()) {
		DA_points >> x >> y >> z;
		Point_3 temp_point = Point_3(x, y, z);
		sampling_points_in_D_A.push_back(temp_point);
		DA_points >> nx >> ny >> nz;
		Point_3 temp_normal = Point_3(nx, ny, nz);
		normal_D_A.push_back(temp_normal);
	}

	ofstream ofile(rootPath + "output\\"+mesh_target+"\\D_S_sampling_points.obj");
	for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
		ofile << "v " << sampling_points_in_D_S[i].x() << " " << sampling_points_in_D_S[i].y() << " " << sampling_points_in_D_S[i].z() << std::endl;
	}
	ofstream ofile_2(rootPath + "output\\"+mesh_target+"\\D_A_sampling_points.obj");
	for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
		ofile_2 << "v " << sampling_points_in_D_A[i].x() << " " << sampling_points_in_D_A[i].y() << " " << sampling_points_in_D_A[i].z() << std::endl;
	}
}

void ReFab::OrientationSamplePoints()
{
	//sample_points_sphere.clear();
	//int number_points = 200;
	//sample_points_sphere.push_back(Point_3(0, 0, 1));
	//for (int i = 1; i <= number_points; i++) {
	//	double phi = acos(-1.0 + (2.0 * i - 1.0) / number_points);
	//	double theta = sqrt(number_points * PI) * phi;
	//	double x = cos(theta) * sin(phi);
	//	double y = sin(theta) * sin(phi);
	//	double z = cos(phi);
	//	if (z >= 0)   //限制采样方向 
	//		sample_points_sphere.push_back(Point_3(x, y, z));
	//}


	//新版高斯球面采样（黄金角度）
	sample_points_sphere.clear();
	int number_points = 200;
	sample_points_sphere.push_back(Point_3(0, 0, 1));

	double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
	double offset = 2.0 / number_points;

	for (int i = 0; i < number_points; ++i) {
		double y = i * offset - 1.0 + (offset / 2.0);
		double radiusXY = sqrt(1.0 - y * y);
		double theta = goldenAngle * i;

		double x = cos(theta) * radiusXY;
		double z = sin(theta) * radiusXY;

		if(z>-0.05)
			sample_points_sphere.push_back({ x * 1, y * 1, z * 1 });
	}
}

std::vector<Point_3> ReFab::OrientationSamplePoints_normal()
{
	vector<Point_3> sample_point_normal;
	sample_point_normal.clear();
	int number_points = 99;  //50 //99
	double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
	double offset = 2.0 / number_points;
	sample_point_normal.push_back({ 0, 0, 1 });
	for (int i = 0; i < number_points; ++i) {
		double y = i * offset - 1.0 + (offset / 2.0);
		double radiusXY = sqrt(1.0 - y * y);
		double theta = goldenAngle * i;

		double x = cos(theta) * radiusXY;
		double z = sin(theta) * radiusXY;
		if (z > 0)
			sample_point_normal.push_back({ x * 1, y * 1, z * 1 });
	}
	return sample_point_normal;
}

std::vector<Point_3> ReFab::OrientationSamplePoints_normal_2()
{
	vector<Point_3> sample_point_normal;
	sample_point_normal.clear();
	int number_points = 99;  //99 //199
	double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
	double offset = 2.0 / number_points;
	sample_point_normal.push_back(Point_3(0, 0, 1));
	for (int i = 0; i < number_points; ++i) {
		double y = i * offset - 1.0 + (offset / 2.0);
		double radiusXY = sqrt(1.0 - y * y);
		double theta = goldenAngle * i;

		double x = cos(theta) * radiusXY;
		double z = sin(theta) * radiusXY;
		sample_point_normal.push_back({ x * 1, y * 1, z * 1 });
	}
	return sample_point_normal;
}

double ReFab::Calculate_area_need_support(Eigen::MatrixXd V_1, Eigen::MatrixXi F_1, double& sum_area)
{
	//ofstream ofile(rootPath + "output\\"+mesh_target+"\\aa.obj");

	double hang_in_area = 0;
	vector<Vector3> vertices_in_face;

	double base_z = 100000;
	for (int i = 0; i < V_1.rows(); i++) {
		//ofile << "v " << V_1(i,0)<< " " << V_1(i, 1) << " " << V_1(i, 2) << std::endl;
		if (V_1(i, 2) < base_z) base_z = V_1(i, 2);
	}
		
	for (int i = 0; i < F_1.rows();i++) {
		vertices_in_face.clear();
		for (int j = 0; j < 3;j++) {
			vertices_in_face.push_back(Vector3(V_1(F_1(i,j),0), V_1(F_1(i, j), 1), V_1(F_1(i, j), 2)));
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
		if((Dot(face_normal, base_normal) + sin(PI / 3.6) >= 0) || (abs(v1.m_z-base_z) < 0.0001))
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
			//ofile << "f " << F_1(i, 0) + 1 << " " << F_1(i, 1) + 1 << " " << F_1(i, 2) + 1 << std::endl;
		}
		sum_area += vec1.cross(vec2).norm() / 2;
	}

	/*for (int i = 0; i < V_1.rows(); i++) {
		if(V_1(i,2) > max_z_of_temp_DI) max_z_of_temp_DI = V_1(i,2);
	}*/

	return hang_in_area;
}

void ReFab::ReOrientation(std::string& meshname1, std::string& meshname2, bool type)
{
	
	vector<Point_3> all_sampling_orientation = OrientationSamplePoints_normal_2();
	/*ofstream ofile2(rootPath + "output\\"+mesh_target+"\\aaaaa.obj");
	for (int i = 0; i < all_sampling_orientation.size(); i++) {
		ofile2 << "v " << all_sampling_orientation[i].x() << " " << all_sampling_orientation[i].y() << " " << all_sampling_orientation[i].z() << std::endl;
	}*/

	Eigen::MatrixXd V_1, V_2;
	Eigen::MatrixXi F_1, F_2;
	igl::readOFF(rootPath + meshname1, V_1, F_1);
	igl::readOFF(rootPath + meshname2, V_2, F_2);
	vector<double> all_volume, all_hang_in_area, all_complexity, all_convex;
	Vector3 center_of_mesh1(0, 0, 0), center_of_mesh2(0, 0, 0);
	double min_z_of_mesh1 = 1000000, min_z_of_mesh2;
	for (int i = 0; i < V_1.rows(); i++) {
		center_of_mesh1.m_x += V_1(i, 0);
		center_of_mesh1.m_y += V_1(i, 1);
		center_of_mesh1.m_z += V_1(i, 2);
		if (V_1(i, 2) < min_z_of_mesh1) min_z_of_mesh1 = V_1(i, 2);
	}
	Vec3 base_point(0, 0, 0);
	int cont_base_point = 0;
	double min_x = 100000, min_y = 100000, max_x = -100000, max_y = -100000;
	for (int i = 0; i < V_1.rows(); i++) {
		if (abs(V_1(i, 2) - min_z_of_mesh1) < 0.0001) {
			base_point.m_x += V_1(i, 0);
			base_point.m_y += V_1(i, 1);
			base_point.m_z += V_1(i, 2);
			if (V_1(i, 0) < min_x) min_x = V_1(i, 0);
			if (V_1(i, 1) < min_y) min_y = V_1(i, 1);
			if (V_1(i, 0) > max_x) max_x = V_1(i, 0);
			if (V_1(i, 1) > max_y) max_y = V_1(i, 1);
			cont_base_point++;
		}
	}
	base_point /= cont_base_point;

	center_of_mesh1 /= V_1.rows();
	for (int i = 0; i < V_2.rows(); i++) {
		center_of_mesh2.m_x += V_2(i, 0);
		center_of_mesh2.m_y += V_2(i, 1);
		center_of_mesh2.m_z += V_2(i, 2);
	}
	center_of_mesh2 /= V_2.rows();

	//旋转，并与mesh1中心对齐
	for (int i = 0; i < V_2.rows(); i++) {
		V_2(i, 0) = V_2(i, 0) - center_of_mesh2.m_x;
		V_2(i, 1) = V_2(i, 1) - center_of_mesh2.m_y;
		V_2(i, 2) = V_2(i, 2) - center_of_mesh2.m_z;
	}
	for (int ori = 0; ori < all_sampling_orientation.size(); ori++) {
		min_z_of_mesh2 = 1000000;
		Eigen::MatrixXd temp_V_2 = V_2;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vector_after(all_sampling_orientation[ori].x(), all_sampling_orientation[ori].y(), all_sampling_orientation[ori].z());
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
		for (int j = 0; j < temp_V_2.rows(); j++) {
			Eigen::MatrixXd temp_V;
			temp_V.resize(3, 1);
			temp_V(0, 0) = temp_V_2(j, 0);
			temp_V(1, 0) = temp_V_2(j, 1);
			temp_V(2, 0) = temp_V_2(j, 2);
			temp_V = rotMatrix.inverse() * temp_V;
			temp_V_2(j, 0) = temp_V(0, 0);
			temp_V_2(j, 1) = temp_V(1, 0);
			temp_V_2(j, 2) = temp_V(2, 0);
		}

		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + center_of_mesh1.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + center_of_mesh1.m_y;
			temp_V_2(i, 2) = temp_V_2(i, 2) + center_of_mesh1.m_z;
		}
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (temp_V_2(i, 2) < min_z_of_mesh2) min_z_of_mesh2 = temp_V_2(i, 2);
		}
		//基底中心对齐
		Vec3 base_point_2(0, 0, 0);
		int cont_base_point_2 = 0;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (abs(temp_V_2(i, 2) - min_z_of_mesh2) < 0.0001) {
				base_point_2.m_x += temp_V_2(i, 0);
				base_point_2.m_y += temp_V_2(i, 1);
				base_point_2.m_z += temp_V_2(i, 2);
				cont_base_point_2++;
			}
		}
		base_point_2 /= cont_base_point_2;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + base_point.m_x - base_point_2.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + base_point.m_y - base_point_2.m_y;
		}
		//基底对齐
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 2) = temp_V_2(i, 2) - min_z_of_mesh2 + min_z_of_mesh1;

			temp_V_2(i, 2) = temp_V_2(i, 2) - 0.01;  //微小偏移，用于防止布尔运算出错  //-0.01
		}

		igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", temp_V_2, F_2);
		Mesh mesh1, mesh2, inter;
		mesh1.clear(); mesh2.clear(); inter.clear();
		CGAL::IO::read_OFF(rootPath + meshname1, mesh1);
		CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", mesh2);
		//CGAL计算相交体积大小
		if (type == 0) {
			Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
			Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
			Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
			Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
			Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
			Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
			if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
				params::vertex_point_map(mesh1_vpm),
				params::vertex_point_map(mesh2_vpm),
				params::vertex_point_map(inter_vpm)))
			{
				std::cout << "Intersection were successfully computed\n";
				CGAL::IO::write_polygon_mesh(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter, CGAL::parameters::stream_precision(17));
			}
			else {
				std::cout << "Intersection could not be computed\n";
				return;
			}
			double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
			all_volume.push_back(volume_temp_DI);
		}
		else {
			Eigen::MatrixXd VA, VB, VC;
			Eigen::VectorXi J, I;
			Eigen::MatrixXi FA, FB, FC;
			igl::MeshBooleanType boolean_type(
				igl::MESH_BOOLEAN_TYPE_INTERSECT);
			igl::readOFF(rootPath + meshname1, VA, FA);
			igl::readOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", VB, FB);
			igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
			igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", VC, FC);
			Mesh inter;
			CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter);
			double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
			all_volume.push_back(volume_temp_DI);
		}

		//计算悬垂区域大小
		double sum_area = 0;
		double hang_in_area = Calculate_area_need_support(temp_V_2, F_2, sum_area);
		all_hang_in_area.push_back(hang_in_area);
		cout << all_volume[all_volume.size() - 1] << " " << all_hang_in_area[all_hang_in_area.size() - 1] << endl;
	}
	double max_volume = -1000000, min_volume = 1000000;
	double max_support = -1000000, min_support = 1000000;
	double max_complexity = -1000000, min_complexity = 1000000;
	double max_convex = -1000000, min_convex = 1000000;
	for (int i = 0; i < all_volume.size(); i++) {
		if (all_hang_in_area[i] > max_support) max_support = all_hang_in_area[i];
		if (all_hang_in_area[i] < min_support) min_support = all_hang_in_area[i];
		if (all_volume[i] > max_volume) max_volume = all_volume[i];
		if (all_volume[i] < min_volume) min_volume = all_volume[i];
	}

	//排序
	for (int i = 0; i < all_hang_in_area.size(); i++)
		for (int j = i + 1; j < all_hang_in_area.size(); j++) {
			double score_i = 0.4 * (all_volume[i] - min_volume) / (max_volume - min_volume) + 0.6 * (1 - (all_hang_in_area[i] - min_support) / (max_support - min_support));
			double score_j = 0.4 * (all_volume[j] - min_volume) / (max_volume - min_volume) + 0.6 * (1 - (all_hang_in_area[j] - min_support) / (max_support - min_support));
			if (score_i < score_j) {
				swap(all_hang_in_area[i], all_hang_in_area[j]);
				swap(all_sampling_orientation[i], all_sampling_orientation[j]);
				swap(all_volume[i], all_volume[j]);
			}
		}


	//*******************************
	//再对最佳方向的几个方向进行平移微调
	double offset_x = (max_x - min_x)/10;
	double offset_y = (max_y - min_y)/10;
	//double offset_x = 2.0;
	//double offset_y = 2.0;

	//double offset = 2.0;  //2.0
	int select_ori = 5;   //选择几个方向进行微调
	int cont_offset = 0;
	vector<double> all_volume_2, all_hang_in_area_2, all_complexity_2, all_convex_2, save_lowest_ratio;
	vector<Point2d> all_offset;
	vector<Point_3> all_sampling_orientation_2;
	for (int ori = 0; ori < select_ori; ori++) {
		min_z_of_mesh2 = 1000000;
		Eigen::MatrixXd temp_V_2 = V_2;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vector_after(all_sampling_orientation[ori].x(), all_sampling_orientation[ori].y(), all_sampling_orientation[ori].z());
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
		for (int j = 0; j < temp_V_2.rows(); j++) {
			Eigen::MatrixXd temp_V;
			temp_V.resize(3, 1);
			temp_V(0, 0) = temp_V_2(j, 0);
			temp_V(1, 0) = temp_V_2(j, 1);
			temp_V(2, 0) = temp_V_2(j, 2);
			temp_V = rotMatrix.inverse() * temp_V;
			temp_V_2(j, 0) = temp_V(0, 0);
			temp_V_2(j, 1) = temp_V(1, 0);
			temp_V_2(j, 2) = temp_V(2, 0);
		}
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + center_of_mesh1.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + center_of_mesh1.m_y;
			temp_V_2(i, 2) = temp_V_2(i, 2) + center_of_mesh1.m_z;
		}
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (temp_V_2(i, 2) < min_z_of_mesh2) min_z_of_mesh2 = temp_V_2(i, 2);
		}

		//基底中心对齐
		Vec3 base_point_2(0, 0, 0);
		int cont_base_point_2 = 0;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (abs(temp_V_2(i, 2) - min_z_of_mesh2) < 0.0001) {
				base_point_2.m_x += temp_V_2(i, 0);
				base_point_2.m_y += temp_V_2(i, 1);
				base_point_2.m_z += temp_V_2(i, 2);
				cont_base_point_2++;
			}
		}
		base_point_2 /= cont_base_point_2;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + base_point.m_x - base_point_2.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + base_point.m_y - base_point_2.m_y;
		}
		//基底对齐
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 2) = temp_V_2(i, 2) - min_z_of_mesh2 + min_z_of_mesh1;
		}

		for (double xx = -offset_x * 2.0; xx <= offset_x * 2.0; xx += offset_x) {
			for (double yy = -offset_y * 2.0; yy <= offset_y * 2.0; yy += offset_y) {
				all_offset.push_back(Point2d(xx, yy));
				cont_offset++;

				Eigen::MatrixXd temp_V_3 = temp_V_2;
				for (int i = 0; i < temp_V_3.rows(); i++) {
					temp_V_3(i, 0) = temp_V_3(i, 0) + xx;
					temp_V_3(i, 1) = temp_V_3(i, 1) + yy;
				}

				igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", temp_V_3, F_2);
				Mesh mesh1, mesh2;
				mesh1.clear(); mesh2.clear();
				CGAL::IO::read_OFF(rootPath + meshname1, mesh1);
				CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", mesh2);
				Mesh inter;
				inter.clear();
				if (type == 0) {
					//CGAL计算相交体积大小
					Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
					Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
					Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
					if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
						params::vertex_point_map(mesh1_vpm),
						params::vertex_point_map(mesh2_vpm),
						params::vertex_point_map(inter_vpm)))
					{
						std::cout << "Intersection were successfully computed\n";
						CGAL::IO::write_polygon_mesh(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter, CGAL::parameters::stream_precision(17));
					}
					else {
						std::cout << "Intersection could not be computed\n";
						return;
					}
					double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
					all_volume_2.push_back(volume_temp_DI);
				}
				else {
					//libigl版本的bool求交，速度慢
					Eigen::MatrixXd VA, VB, VC;
					Eigen::VectorXi J, I;
					Eigen::MatrixXi FA, FB, FC;
					igl::MeshBooleanType boolean_type(
						igl::MESH_BOOLEAN_TYPE_INTERSECT);
					igl::readOFF(rootPath + meshname1, VA, FA);
					igl::readOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", VB, FB);
					igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
					igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", VC, FC);
					CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter);
					double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
					all_volume_2.push_back(volume_temp_DI);
				}


				all_sampling_orientation_2.push_back(all_sampling_orientation[ori]);
				//all_volume_2.push_back(volume_temp_DI);
				cout << "orientation: " << cont_offset - 1 << endl;
				//计算悬垂区域大小
				double sum_area = 0;
				double max_z_of_temp_DI = -100000;
				double hang_in_area = Calculate_area_need_support(temp_V_3, F_2, sum_area);
				for (int i = 0; i < V_1.rows(); i++)
					if (V_1(i, 2) > max_z_of_temp_DI) max_z_of_temp_DI = V_1(i, 2);

				double mean_mean_principal_curvature = 0, abs_mean_mean_principal_curvature = 0;
				double sum_additive_interface_area = 0;
				double min_z_of_additive_interface = 100000;
				//这个函数需要加速！！！(占了60％的时间)
				Get_additive_inter_face(inter, mesh2, mesh1, rootPath + "output\\" + mesh_target + "\\temp_Interface_A.off",
					sum_additive_interface_area, min_z_of_additive_interface);
				//Calculate_mean_curvature(rootPath + "output\\" + mesh_target + "\\temp_Interface_A.off", false, mean_mean_principal_curvature, abs_mean_mean_principal_curvature);
				double lowest_position_ratio = (min_z_of_additive_interface - min_z_of_mesh1) / (max_z_of_temp_DI - min_z_of_mesh1);
				//all_convex_2.push_back(sqrt(mean_mean_principal_curvature) * lowest_position_ratio * sqrt(1 - (sum_additive_interface_area / sum_area)));
// 
				all_convex_2.push_back(lowest_position_ratio + pow((1 - (sum_additive_interface_area / sum_area)), 2));
				save_lowest_ratio.push_back(lowest_position_ratio);
				//all_convex_2.push_back(lowest_position_ratio);

				//all_convex_2.push_back(lowest_position_ratio);
				cout << all_volume_2[all_volume_2.size() - 1] << " " << lowest_position_ratio << endl;
			}
		}
	}

	max_volume = -1000000, min_volume = 1000000;
	max_support = -1000000, min_support = 1000000;
	max_complexity = -1000000, min_complexity = 1000000;
	max_convex = -1000000, min_convex = 1000000;
	for (int i = 0; i < all_volume_2.size(); i++) {
		if (all_volume_2[i] > max_volume) max_volume = all_volume_2[i];
		if (all_volume_2[i] < min_volume) min_volume = all_volume_2[i];
		if (all_convex_2[i] > max_convex) max_convex = all_convex_2[i];
		if (all_convex_2[i] < min_convex) min_convex = all_convex_2[i];
	}

	double max_score_2 = -1000000;
	int max_index_2 = 0;
	for (int i = 0; i < all_volume_2.size(); i++) {
		double score_volume = (all_volume_2[i] - min_volume) / (max_volume - min_volume);
		//double score_complexity = 1 - (all_complexity_2[i] - min_complexity) / (max_complexity - min_complexity);
		double score_complexity = 0;
		double score_convex = (all_convex_2[i] - min_convex) / (max_convex - min_convex);
		double score = 0.4 * score_volume + 0.6 * score_convex;
		if (score > max_score_2) {
			max_score_2 = score;
			max_index_2 = i;
		}
	}
	record_data[1] = save_lowest_ratio[max_index_2];

	//重跑一次
	cout << "Best orientation: " << max_index_2 << endl;
	Eigen::MatrixXd temp_V_2 = V_2;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(all_sampling_orientation_2[max_index_2].x(), all_sampling_orientation_2[max_index_2].y(), all_sampling_orientation_2[max_index_2].z());
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int j = 0; j < temp_V_2.rows(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = temp_V_2(j, 0);
		temp_V(1, 0) = temp_V_2(j, 1);
		temp_V(2, 0) = temp_V_2(j, 2);
		temp_V = rotMatrix.inverse() * temp_V;
		temp_V_2(j, 0) = temp_V(0, 0);
		temp_V_2(j, 1) = temp_V(1, 0);
		temp_V_2(j, 2) = temp_V(2, 0);
	}
	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 0) = temp_V_2(i, 0) + center_of_mesh1.m_x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + center_of_mesh1.m_y;
		temp_V_2(i, 2) = temp_V_2(i, 2) + center_of_mesh1.m_z;
	}
	min_z_of_mesh2 = 1000000;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		if (temp_V_2(i, 2) < min_z_of_mesh2) min_z_of_mesh2 = temp_V_2(i, 2);
	}

	//基底中心对齐
	Vec3 base_point_2(0, 0, 0);
	int cont_base_point_2 = 0;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		if (abs(temp_V_2(i, 2) - min_z_of_mesh2) < 0.0001) {
			base_point_2.m_x += temp_V_2(i, 0);
			base_point_2.m_y += temp_V_2(i, 1);
			base_point_2.m_z += temp_V_2(i, 2);
			cont_base_point_2++;
		}
	}
	base_point_2 /= cont_base_point_2;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 0) = temp_V_2(i, 0) + base_point.m_x - base_point_2.m_x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + base_point.m_y - base_point_2.m_y;
	}

	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 2) = temp_V_2(i, 2) - min_z_of_mesh2 + min_z_of_mesh1;
		temp_V_2(i, 0) = temp_V_2(i, 0) + all_offset[max_index_2].x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + all_offset[max_index_2].y;
		temp_V_2(i, 0) = temp_V_2(i, 0);
		temp_V_2(i, 1) = temp_V_2(i, 1);
	}
	igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", temp_V_2, F_2);
	Mesh mesh1, mesh2, inter;
	mesh1.clear(); mesh2.clear(); inter.clear();
	CGAL::IO::read_OFF(rootPath + meshname1, mesh1);
	CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", mesh2);
	Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
	Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
	Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
	if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
		params::vertex_point_map(mesh1_vpm),
		params::vertex_point_map(mesh2_vpm),
		params::vertex_point_map(inter_vpm)))
	{
		std::cout << "Intersection were successfully computed\n";
		CGAL::IO::write_polygon_mesh(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter, CGAL::parameters::stream_precision(17));
	}
	igl::writeOFF(rootPath + meshname2.substr(0, meshname2.size() - 4) + "_new.off", temp_V_2, F_2);
	std::cout << "final ori:" << all_sampling_orientation_2[max_index_2].x() << " " << all_sampling_orientation_2[max_index_2].y() << " " << all_sampling_orientation_2[max_index_2].z() << endl;
	std::cout << "final offset:" << all_offset[max_index_2].x << " " << all_offset[max_index_2].y << endl;
}

void ReFab::ReOrientation_combined_optimization(std::string& meshname1, std::string& meshname2, bool type)
{

	vector<Point_3> all_sampling_orientation = OrientationSamplePoints_normal_2();
	/*ofstream ofile2(rootPath + "output\\"+mesh_target+"\\aaaaa.obj");
	for (int i = 0; i < all_sampling_orientation.size(); i++) {
		ofile2 << "v " << all_sampling_orientation[i].x() << " " << all_sampling_orientation[i].y() << " " << all_sampling_orientation[i].z() << std::endl;
	}*/

	Eigen::MatrixXd V_1, V_2;
	Eigen::MatrixXi F_1, F_2;
	igl::readOFF(rootPath + meshname1, V_1, F_1);
	igl::readOFF(rootPath + meshname2, V_2, F_2);
	vector<double> all_volume, all_hang_in_area, all_complexity, all_convex;
	Vector3 center_of_mesh1(0, 0, 0), center_of_mesh2(0, 0, 0);
	double min_z_of_mesh1 = 1000000, min_z_of_mesh2;
	for (int i = 0; i < V_1.rows(); i++) {
		center_of_mesh1.m_x += V_1(i, 0);
		center_of_mesh1.m_y += V_1(i, 1);
		center_of_mesh1.m_z += V_1(i, 2);
		if (V_1(i, 2) < min_z_of_mesh1) min_z_of_mesh1 = V_1(i, 2);
	}
	Vec3 base_point(0, 0, 0);
	int cont_base_point = 0;
	double min_x = 100000, min_y = 100000, max_x = -100000, max_y = -100000;
	for (int i = 0; i < V_1.rows(); i++) {
		if (abs(V_1(i, 2) - min_z_of_mesh1) < 0.0001) {
			base_point.m_x += V_1(i, 0);
			base_point.m_y += V_1(i, 1);
			base_point.m_z += V_1(i, 2);
			if (V_1(i, 0) < min_x) min_x = V_1(i, 0);
			if (V_1(i, 1) < min_y) min_y = V_1(i, 1);
			if (V_1(i, 0) > max_x) max_x = V_1(i, 0);
			if (V_1(i, 1) > max_y) max_y = V_1(i, 1);
			cont_base_point++;
		}
	}
	base_point /= cont_base_point;

	center_of_mesh1 /= V_1.rows();
	for (int i = 0; i < V_2.rows(); i++) {
		center_of_mesh2.m_x += V_2(i, 0);
		center_of_mesh2.m_y += V_2(i, 1);
		center_of_mesh2.m_z += V_2(i, 2);
	}
	center_of_mesh2 /= V_2.rows();

	//旋转，并与mesh1中心对齐
	for (int i = 0; i < V_2.rows(); i++) {
		V_2(i, 0) = V_2(i, 0) - center_of_mesh2.m_x;
		V_2(i, 1) = V_2(i, 1) - center_of_mesh2.m_y;
		V_2(i, 2) = V_2(i, 2) - center_of_mesh2.m_z;
	}




	//*******************************
	//再对最佳方向的几个方向进行平移微调
	double offset_x = (max_x - min_x) / 10;
	double offset_y = (max_y - min_y) / 10; 
	if (offset_x == 0)
		offset_x = 1;
	if (offset_y == 0)
		offset_y = 1;
	//double offset_x = 2.0;
	//double offset_y = 2.0;

	//double offset = 2.0;  //2.0
	int select_ori = all_sampling_orientation.size();   //选择几个方向进行微调
	int cont_offset = 0;
	vector<double> all_volume_2, all_hang_in_area_2, all_complexity_2, all_convex_2, save_lowest_ratio;
	vector<Point2d> all_offset;
	vector<Point_3> all_sampling_orientation_2;
	for (int ori = 0; ori < select_ori; ori++) {
		min_z_of_mesh2 = 1000000;
		Eigen::MatrixXd temp_V_2 = V_2;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vector_after(all_sampling_orientation[ori].x(), all_sampling_orientation[ori].y(), all_sampling_orientation[ori].z());
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
		for (int j = 0; j < temp_V_2.rows(); j++) {
			Eigen::MatrixXd temp_V;
			temp_V.resize(3, 1);
			temp_V(0, 0) = temp_V_2(j, 0);
			temp_V(1, 0) = temp_V_2(j, 1);
			temp_V(2, 0) = temp_V_2(j, 2);
			temp_V = rotMatrix.inverse() * temp_V;
			temp_V_2(j, 0) = temp_V(0, 0);
			temp_V_2(j, 1) = temp_V(1, 0);
			temp_V_2(j, 2) = temp_V(2, 0);
		}
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + center_of_mesh1.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + center_of_mesh1.m_y;
			temp_V_2(i, 2) = temp_V_2(i, 2) + center_of_mesh1.m_z;
		}
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (temp_V_2(i, 2) < min_z_of_mesh2) min_z_of_mesh2 = temp_V_2(i, 2);
		}

		//基底中心对齐
		Vec3 base_point_2(0, 0, 0);
		int cont_base_point_2 = 0;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			if (abs(temp_V_2(i, 2) - min_z_of_mesh2) < 0.0001) {
				base_point_2.m_x += temp_V_2(i, 0);
				base_point_2.m_y += temp_V_2(i, 1);
				base_point_2.m_z += temp_V_2(i, 2);
				cont_base_point_2++;
			}
		}
		base_point_2 /= cont_base_point_2;
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 0) = temp_V_2(i, 0) + base_point.m_x - base_point_2.m_x;
			temp_V_2(i, 1) = temp_V_2(i, 1) + base_point.m_y - base_point_2.m_y;
		}
		//基底对齐
		for (int i = 0; i < temp_V_2.rows(); i++) {
			temp_V_2(i, 2) = temp_V_2(i, 2) - min_z_of_mesh2 + min_z_of_mesh1;
		}

		for (double xx = -offset_x * 2.0; xx <= offset_x * 2.0; xx += offset_x/2) {  //xx += offset_x
			for (double yy = -offset_y * 2.0; yy <= offset_y * 2.0; yy += offset_y/2) {
				all_offset.push_back(Point2d(xx, yy));
				cont_offset++;
				//if (ori == 0 || ori==90 || ori==91)
				//	continue;

				Eigen::MatrixXd temp_V_3 = temp_V_2;
				for (int i = 0; i < temp_V_3.rows(); i++) {
					temp_V_3(i, 0) = temp_V_3(i, 0) + xx;
					temp_V_3(i, 1) = temp_V_3(i, 1) + yy;
				}

				igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", temp_V_3, F_2);
				Mesh mesh1, mesh2;
				mesh1.clear(); mesh2.clear();
				CGAL::IO::read_OFF(rootPath + meshname1, mesh1);
				CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", mesh2);
				Mesh inter;
				inter.clear();
				if (type == 0) {
					//CGAL计算相交体积大小
					Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
					Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
					Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
					Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
					if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
						params::vertex_point_map(mesh1_vpm),
						params::vertex_point_map(mesh2_vpm),
						params::vertex_point_map(inter_vpm)))
					{
						std::cout << "Intersection were successfully computed\n";
						CGAL::IO::write_polygon_mesh(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter, CGAL::parameters::stream_precision(17));
					}
					else {
						std::cout << "Intersection could not be computed\n";
						return;
					}
					double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
					all_volume_2.push_back(volume_temp_DI);
				}
				else {
					//libigl版本的bool求交，速度慢
					Eigen::MatrixXd VA, VB, VC;
					Eigen::VectorXi J, I;
					Eigen::MatrixXi FA, FB, FC;
					igl::MeshBooleanType boolean_type(
						igl::MESH_BOOLEAN_TYPE_INTERSECT);
					igl::readOFF(rootPath + meshname1, VA, FA);
					igl::readOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", VB, FB);
					igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
					igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", VC, FC);
					CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter);
					double volume_temp_DI = CGAL::Polygon_mesh_processing::volume(inter);
					all_volume_2.push_back(volume_temp_DI);
				}


				all_sampling_orientation_2.push_back(all_sampling_orientation[ori]);
				//all_volume_2.push_back(volume_temp_DI);
				cout << "orientation: " << ori << endl;
				cout << "number: " << cont_offset - 1 << endl;
				//计算悬垂区域大小
				double sum_area = 0;
				double max_z_of_temp_DI = -100000;
				double hang_in_area = Calculate_area_need_support(temp_V_3, F_2, sum_area);
				for (int i = 0; i < V_1.rows(); i++)
					if (V_1(i, 2) > max_z_of_temp_DI) max_z_of_temp_DI = V_1(i, 2);

				double mean_mean_principal_curvature = 0, abs_mean_mean_principal_curvature = 0;
				double sum_additive_interface_area = 0;
				double min_z_of_additive_interface = 100000;
				//这个函数需要加速！！！(占了60％的时间)
				Get_additive_inter_face(inter, mesh2, mesh1, rootPath + "output\\" + mesh_target + "\\temp_Interface_A.off",
					sum_additive_interface_area, min_z_of_additive_interface);
				//Calculate_mean_curvature(rootPath + "output\\" + mesh_target + "\\temp_Interface_A.off", false, mean_mean_principal_curvature, abs_mean_mean_principal_curvature);
				double lowest_position_ratio = (min_z_of_additive_interface - min_z_of_mesh1) / (max_z_of_temp_DI - min_z_of_mesh1);
				//all_convex_2.push_back(sqrt(mean_mean_principal_curvature) * lowest_position_ratio * sqrt(1 - (sum_additive_interface_area / sum_area)));
// 
				all_convex_2.push_back(lowest_position_ratio + pow((1 - (sum_additive_interface_area / sum_area)), 2));
				save_lowest_ratio.push_back(lowest_position_ratio);
				//all_convex_2.push_back(lowest_position_ratio);

				//all_convex_2.push_back(lowest_position_ratio);
				cout << all_volume_2[all_volume_2.size() - 1] << " " << lowest_position_ratio << endl;
			}
		}
	}

	double max_volume = -1000000, min_volume = 1000000;
	double max_support = -1000000, min_support = 1000000;
	double max_complexity = -1000000, min_complexity = 1000000;
	double max_convex = -1000000, min_convex = 1000000;
	for (int i = 0; i < all_volume_2.size(); i++) {
		if (all_volume_2[i] > max_volume) max_volume = all_volume_2[i];
		if (all_volume_2[i] < min_volume) min_volume = all_volume_2[i];
		if (all_convex_2[i] > max_convex) max_convex = all_convex_2[i];
		if (all_convex_2[i] < min_convex) min_convex = all_convex_2[i];
	}

	double max_score_2 = -1000000;
	int max_index_2 = 0;
	for (int i = 0; i < all_volume_2.size(); i++) {
		double score_volume = (all_volume_2[i] - min_volume) / (max_volume - min_volume);
		//double score_complexity = 1 - (all_complexity_2[i] - min_complexity) / (max_complexity - min_complexity);
		double score_complexity = 0;
		double score_convex = (all_convex_2[i] - min_convex) / (max_convex - min_convex);
		double score = 0.4 * score_volume + 0.6 * score_convex;
		if (score > max_score_2) {
			max_score_2 = score;
			max_index_2 = i;
		}
	}
	record_data[1] = save_lowest_ratio[max_index_2];

	//重跑一次
	cout << "Best orientation: " << max_index_2 << endl;
	Eigen::MatrixXd temp_V_2 = V_2;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(all_sampling_orientation_2[max_index_2].x(), all_sampling_orientation_2[max_index_2].y(), all_sampling_orientation_2[max_index_2].z());
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int j = 0; j < temp_V_2.rows(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = temp_V_2(j, 0);
		temp_V(1, 0) = temp_V_2(j, 1);
		temp_V(2, 0) = temp_V_2(j, 2);
		temp_V = rotMatrix.inverse() * temp_V;
		temp_V_2(j, 0) = temp_V(0, 0);
		temp_V_2(j, 1) = temp_V(1, 0);
		temp_V_2(j, 2) = temp_V(2, 0);
	}
	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 0) = temp_V_2(i, 0) + center_of_mesh1.m_x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + center_of_mesh1.m_y;
		temp_V_2(i, 2) = temp_V_2(i, 2) + center_of_mesh1.m_z;
	}
	min_z_of_mesh2 = 1000000;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		if (temp_V_2(i, 2) < min_z_of_mesh2) min_z_of_mesh2 = temp_V_2(i, 2);
	}

	//基底中心对齐
	Vec3 base_point_2(0, 0, 0);
	int cont_base_point_2 = 0;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		if (abs(temp_V_2(i, 2) - min_z_of_mesh2) < 0.0001) {
			base_point_2.m_x += temp_V_2(i, 0);
			base_point_2.m_y += temp_V_2(i, 1);
			base_point_2.m_z += temp_V_2(i, 2);
			cont_base_point_2++;
		}
	}
	base_point_2 /= cont_base_point_2;
	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 0) = temp_V_2(i, 0) + base_point.m_x - base_point_2.m_x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + base_point.m_y - base_point_2.m_y;
	}

	for (int i = 0; i < temp_V_2.rows(); i++) {
		temp_V_2(i, 2) = temp_V_2(i, 2) - min_z_of_mesh2 + min_z_of_mesh1;
		temp_V_2(i, 0) = temp_V_2(i, 0) + all_offset[max_index_2].x;
		temp_V_2(i, 1) = temp_V_2(i, 1) + all_offset[max_index_2].y;
		temp_V_2(i, 0) = temp_V_2(i, 0);
		temp_V_2(i, 1) = temp_V_2(i, 1);
	}
	igl::writeOFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", temp_V_2, F_2);
	Mesh mesh1, mesh2, inter;
	mesh1.clear(); mesh2.clear(); inter.clear();
	CGAL::IO::read_OFF(rootPath + meshname1, mesh1);
	CGAL::IO::read_OFF(rootPath + "output\\" + mesh_target + "\\temp_target.off", mesh2);
	Exact_point_map mesh1_exact_points = mesh1.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_point_map mesh2_exact_points = mesh2.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_point_map inter_exact_points = inter.add_property_map<vertex_descriptor, EK::Point_3>("v:exact_point").first;
	Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
	Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);
	Exact_vertex_point_map inter_vpm(inter_exact_points, inter);
	if (PMP::corefine_and_compute_intersection(mesh1, mesh2, inter,
		params::vertex_point_map(mesh1_vpm),
		params::vertex_point_map(mesh2_vpm),
		params::vertex_point_map(inter_vpm)))
	{
		std::cout << "Intersection were successfully computed\n";
		CGAL::IO::write_polygon_mesh(rootPath + "output\\" + mesh_target + "\\temp_DI.off", inter, CGAL::parameters::stream_precision(17));
	}
	igl::writeOFF(rootPath + meshname2.substr(0, meshname2.size() - 4) + "_new.off", temp_V_2, F_2);
	std::cout << "final ori:" << all_sampling_orientation_2[max_index_2].x() << " " << all_sampling_orientation_2[max_index_2].y() << " " << all_sampling_orientation_2[max_index_2].z() << endl;
	std::cout << "final offset:" << all_offset[max_index_2].x << " " << all_offset[max_index_2].y << endl;
}

void ReFab::Get_additive_inter_face(Mesh inter, Mesh mesh_a, Mesh mesh_s, std::string output_a, double& sum_additive_interface_area, double& min_z_of_additive_interface)
{
	//如果一个face的中点在mesh_a内部，那么这个点是增材接触面，不考虑减材接触面
	Mesh interface_a;

	//建立一个新的映射，对于interface_s中的每一个点，都有一个对应的点在inter中
	std::unordered_map<Mesh::Vertex_index, Mesh::Vertex_index> interface_a_map;

	//创建一个内外判别器，用于判断点是否在mesh_a内部
	CGAL::Side_of_triangle_mesh<Mesh, K> inside_a(mesh_a);
	//遍历inter中的每一个面
	for (Mesh::Face_index face : faces(inter))
	{
		//就只计算一个面的终点是不是在mesh_a内部

		double x = 0.0, y = 0.0, z = 0.0;
		std::vector<Point_3> PointOfFace;

		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			Point_3 p = inter.point(vd);
			PointOfFace.push_back(p);
			//将xyz累加起来
			x += p.x();
			y += p.y();
			z += p.z();
		}
		//计算face的法向
		K::Vector_3 normal = CGAL::normal(PointOfFace[0], PointOfFace[1], PointOfFace[2]);
		K::Vector_3 center(x / 3.0, y / 3.0, z / 3.0);
		center = center + 0.01 * normal;


		Point_3 centerP(center.x(), center.y(), center.z());
		//判断中心点是否在mesh_a内部
		if (inside_a(centerP) == CGAL::ON_CONVEX_SIDE) {//说明这个面落在mesh_a内部
			//存储该面对应的点的索引
			std::vector<Mesh::Vertex_index> new_vertices;
			//将组成面的点保存到interface_a中
			for (Mesh::Vertex_index vertex_index : vertices_around_face(inter.halfedge(face), inter)) {
				//如果顶点还没有添加到interface_s，则添加它
				if (interface_a_map.find(vertex_index) == interface_a_map.end()) {
					Point_3 p = inter.point(vertex_index);
					Mesh::Vertex_index new_vertex = interface_a.add_vertex(p);
					interface_a_map[vertex_index] = new_vertex;
				}
				new_vertices.push_back(interface_a_map[vertex_index]);
			}

			//添加新的面到interface_s中,对于定点数不同的面分别处理
			if (new_vertices.size() == 3)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2]);
			else if (new_vertices.size() == 4)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2], new_vertices[3]);
			else
				std::cerr << "Error: Unable to add face to interface_a." << std::endl;

		}
	}
	//将interface_s保存到output1中
	std::cout << "Split finished\n";
	//输出增材接触面
	CGAL::IO::write_polygon_mesh(output_a, interface_a, CGAL::parameters::stream_precision(17));

	//计算增材接触面的面积以及最小z值
	for (Mesh::Face_index face : faces(interface_a)) {
		Eigen::Vector3d vec1, vec2;
		std::vector<Point_3> PointOfFace;
		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			Point_3 p = inter.point(vd);
			PointOfFace.push_back(p);
		}
		Vector3 v1(PointOfFace[0].x(), PointOfFace[0].y(), PointOfFace[0].z());
		Vector3 v2(PointOfFace[1].x(), PointOfFace[1].y(), PointOfFace[1].z());
		Vector3 v3(PointOfFace[2].x(), PointOfFace[2].y(), PointOfFace[2].z());
		vec1.x() = v2.m_x - v1.m_x;
		vec1.y() = v2.m_y - v1.m_y;
		vec1.z() = v2.m_z - v1.m_z;
		vec2.x() = v3.m_x - v1.m_x;
		vec2.y() = v3.m_y - v1.m_y;
		vec2.z() = v3.m_z - v1.m_z;
		sum_additive_interface_area += vec1.cross(vec2).norm() / 2;
	}
	for (Mesh::Vertex_index vd : vertices(interface_a)) {
		if (interface_a.point(vd).z() < min_z_of_additive_interface) min_z_of_additive_interface = interface_a.point(vd).z();
	}
}

void ReFab::Get_additive_inter_face(std::string input_inter, std::string input_a, std::string input_s, std::string output_a, double& sum_additive_interface_area, double& min_z_of_additive_interface)
{

	Mesh inter,mesh_a,mesh_s;
	PMP::IO::read_polygon_mesh(input_inter, inter);
	PMP::IO::read_polygon_mesh(input_a, mesh_a);
	PMP::IO::read_polygon_mesh(input_s, mesh_s);

	Mesh interface_a;
	std::unordered_map<Mesh::Vertex_index, Mesh::Vertex_index> interface_a_map;

	//创建一个内外判别器，用于判断点是否在mesh1内部
	CGAL::Side_of_triangle_mesh<Mesh, K> inside_s(mesh_s), inside_a(mesh_a);
	//遍历inter中的每一个面
	double min_z_of_DI = 100000;
	for (Mesh::Face_index face : faces(inter)) {  
		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			if (inter.point(vd).z() < min_z_of_DI)
				min_z_of_DI = inter.point(vd).z();
		}
	}

	for (Mesh::Face_index face : faces(inter))
	{
		//不考虑底部的面（计算时可能会出现误差）
		bool jud_continue = false;
		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			if (abs(inter.point(vd).z() - min_z_of_DI) < 0.0001) {
				jud_continue = true;
				break;
			}
		}
		if (jud_continue == true) continue;


		//如果inter的一个面的中心点在mesh1内部，而且还必须有一个顶点在mesh1的内部，则将inter保存到interface_s中
		bool is_in_mesh_s = false;

		//对于一个面，它的顶点落在mesh1内部的个数，和落在mesh2内部的个数
		int cnt_ins = 0, cnt_ina = 0, cnt_ons = 0, cnt_ona = 0;
		std::vector<Point_3> PointOfFace;
		for (Mesh::Vertex_index vd : inter.vertices_around_face(inter.halfedge(face))) {
			if (inter.point(vd).z() < min_z_of_DI)
				min_z_of_DI = inter.point(vd).z();
			Point_3 p = inter.point(vd);
			PointOfFace.push_back(p);
			CGAL::Bounded_side res = inside_s(p);
			if (res == CGAL::ON_BOUNDED_SIDE) cnt_ins++;
			if (res == CGAL::ON_BOUNDARY) cnt_ons++;
			CGAL::Bounded_side res2 = inside_a(p);
			if (res == CGAL::ON_BOUNDED_SIDE) cnt_ina++;
			if (res2 == CGAL::ON_BOUNDARY) cnt_ona++;
		}


		//只需要判断在a和s表面上点的个数就可以判断了
		if (cnt_ona >= cnt_ons) is_in_mesh_s = true;

		//这个面在mesh1内部
		//  !is_in_mesh_s && cnt_ona <3 && cnt_ons != 2
		//  !is_in_mesh_s
		//  cnt_ins<2 && cnt_ins<2
		if(!is_in_mesh_s) {
			//存储该面对应的点的索引
			std::vector<Mesh::Vertex_index> new_vertices;
			//将组成面的点保存到interface_a中
			for (Mesh::Vertex_index vertex_index : vertices_around_face(inter.halfedge(face), inter)) {
				//如果顶点还没有添加到interface_s，则添加它
				if (interface_a_map.find(vertex_index) == interface_a_map.end()) {
					Point_3 p = inter.point(vertex_index);
					Mesh::Vertex_index new_vertex = interface_a.add_vertex(p);
					interface_a_map[vertex_index] = new_vertex;
				}
				new_vertices.push_back(interface_a_map[vertex_index]);
			}

			//添加新的面到interface_s中,对于定点数不同的面分别处理
			if (new_vertices.size() == 3)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2]);
			else if (new_vertices.size() == 4)
				Mesh::Face_index new_face = interface_a.add_face(new_vertices[0], new_vertices[1], new_vertices[2], new_vertices[3]);
			else
				std::cerr << "Error: Unable to add face to interface_a." << std::endl;

		}
	}

	//将interface_s保存到output1中
	//std::cout << "Split finished\n";
	CGAL::IO::write_polygon_mesh(output_a, interface_a, CGAL::parameters::stream_precision(17));


	Eigen::MatrixXd V_f;
	Eigen::MatrixXi F_f;
	igl::readOFF(output_a, V_f, F_f);
	for (int i = 0; i < F_f.rows(); i++) {
		Eigen::Vector3d vec1, vec2;
		Vector3 v1(V_f(F_f(i, 0), 0), V_f(F_f(i, 0), 1), V_f(F_f(i, 0), 2));
		Vector3 v2(V_f(F_f(i, 1), 0), V_f(F_f(i, 1), 1), V_f(F_f(i, 1), 2));
		Vector3 v3(V_f(F_f(i, 2), 0), V_f(F_f(i, 2), 1), V_f(F_f(i, 2), 2));
		vec1.x() = v2.m_x - v1.m_x;
		vec1.y() = v2.m_y - v1.m_y;
		vec1.z() = v2.m_z - v1.m_z;
		vec2.x() = v3.m_x - v1.m_x;
		vec2.y() = v3.m_y - v1.m_y;
		vec2.z() = v3.m_z - v1.m_z;
		sum_additive_interface_area += vec1.cross(vec2).norm() / 2;
	}
	for(int i =0;i<V_f.rows();i++){
		if(V_f(i,2) < min_z_of_additive_interface) min_z_of_additive_interface = V_f(i,2);
	}
}

Point_3 ReFab::calculate_normal(Point_3 v1, Point_3 v2, Point_3 v3)
{
	double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
	double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
	double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());
	Vec3 vn(na, nb, nc);
	vn.Normalized();
	Point_3 vn_2(vn.m_x, vn.m_y, vn.m_z);
	return vn_2;
}

void ReFab::Collision_Detection_For_DS_and_DA()
{

	OrientationSamplePoints();
	/*ofstream ttofile(rootPath + "output\\"+mesh_target+"\\bbb.obj");
	for (int i = 0; i < sample_points_sphere.size(); i++)
	{
		ttofile << "v " << sample_points_sphere[i].x() << " " << sample_points_sphere[i].y() << " " << sample_points_sphere[i].z() << std::endl;
	}*/

	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_covering_points_S, flag_covering_points_A;

	flag_accessible_points_S.resize(sampling_points_in_D_S.size());
	flag_accessible_points_A.resize(sampling_points_in_D_A.size());
	flag_covering_points_S.resize(sampling_points_in_D_I.size());
	flag_covering_points_A.resize(sampling_points_in_D_I.size());

	for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
		flag_accessible_points_S[i] = false;
	}
	for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
		flag_accessible_points_A[i] = false;
	}
	for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
		flag_covering_points_S[i] = false;
		flag_covering_points_A[i] = false;
	}

	//碰撞检测，计算不可达点（任意方向都与Di上的点碰撞）
	for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
		///////////////////rotate/////////////////////
		std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I,temp_normal_S,temp_normal_A;
		temp_V_S.resize(sampling_points_in_D_S.size());
		temp_V_A.resize(sampling_points_in_D_A.size());
		temp_V_I.resize(sampling_points_in_D_I.size());
		temp_normal_S.resize(normal_D_S.size());
		temp_normal_A.resize(normal_D_A.size());
		for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
			temp_V_S[i].resize(3, 1);
			temp_V_S[i](0, 0) = sampling_points_in_D_S[i].x();
			temp_V_S[i](1, 0) = sampling_points_in_D_S[i].y();
			temp_V_S[i](2, 0) = sampling_points_in_D_S[i].z();
		}
		for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
			temp_V_A[i].resize(3, 1);
			temp_V_A[i](0, 0) = sampling_points_in_D_A[i].x();
			temp_V_A[i](1, 0) = sampling_points_in_D_A[i].y();
			temp_V_A[i](2, 0) = sampling_points_in_D_A[i].z();
		}
		for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
			temp_V_I[i].resize(3, 1);
			temp_V_I[i](0, 0) = sampling_points_in_D_I[i].x();
			temp_V_I[i](1, 0) = sampling_points_in_D_I[i].y();
			temp_V_I[i](2, 0) = sampling_points_in_D_I[i].z();
		}
		for (int i = 0; i < normal_D_S.size(); i++) {
			temp_normal_S[i].resize(3, 1);
			temp_normal_S[i](0, 0) = normal_D_S[i].x();
			temp_normal_S[i](1, 0) = normal_D_S[i].y();
			temp_normal_S[i](2, 0) = normal_D_S[i].z();
		}
		for (int i = 0; i < normal_D_A.size(); i++) {
			temp_normal_A[i].resize(3, 1);
			temp_normal_A[i](0, 0) = normal_D_A[i].x();
			temp_normal_A[i](1, 0) = normal_D_A[i].y();
			temp_normal_A[i](2, 0) = normal_D_A[i].z();
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < sampling_points_in_D_S.size(); i++)
			temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
		for (int i = 0; i < sampling_points_in_D_A.size(); i++)
			temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
		for (int i = 0; i < sampling_points_in_D_I.size(); i++)
			temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];
		for (int i = 0; i < normal_D_S.size(); i++)
			temp_normal_S[i] = rotMatrix.inverse() * temp_normal_S[i];
		for (int i = 0; i < normal_D_A.size(); i++)
			temp_normal_A[i] = rotMatrix.inverse() * temp_normal_A[i];
		//////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points_A = 0, cont_unaccessible_points_A = 0, cont_accessible_points_S = 0, cont_unaccessible_points_S = 0;

		for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
			if (flag_accessible_points_S[i] == true) {
				cont_accessible_points_S++;
				continue;
			}
			cont_unaccessible_points_S++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_S[i](0, 0) + Cutter.ball_r * temp_normal_S[i](0, 0);
			center_point.m_y = temp_V_S[i](1, 0) + Cutter.ball_r * temp_normal_S[i](1, 0);
			center_point.m_z = temp_V_S[i](2, 0) + Cutter.ball_r * temp_normal_S[i](2, 0);
			Vector3 center_point_2;
			center_point_2.m_x = center_point_2.m_y = center_point_2.m_z = 0;
			center_point_2.m_x = temp_V_S[i](0, 0) - Cutter.ball_r * temp_normal_S[i](0, 0);
			center_point_2.m_y = temp_V_S[i](1, 0) - Cutter.ball_r * temp_normal_S[i](1, 0);
			center_point_2.m_z = temp_V_S[i](2, 0) - Cutter.ball_r * temp_normal_S[i](2, 0);
			for (int ii = 0; ii < sampling_points_in_D_I.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if ((temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r) || (temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r)) {
					continue;
				}
				if ((temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) && (temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height)) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				flag_accessible_points_S[i] = true;
				cont_accessible_points_S++;
				cont_unaccessible_points_S--;
			}
		}
		std::cout << "id of orientation:" << ori << std::endl;
		std::cout << "#accessible points of D_S:" << cont_accessible_points_S << std::endl;
		std::cout << "#unaccessible points of D_S:" << cont_unaccessible_points_S << std::endl << std::endl;

		for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
			if (flag_accessible_points_A[i] == true) {
				cont_accessible_points_A++;
				continue;
			}
			cont_unaccessible_points_A++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_A[i](0, 0);
			center_point.m_y = temp_V_A[i](1, 0);
			center_point.m_z = temp_V_A[i](2, 0);
			for (int ii = 0; ii < sampling_points_in_D_I.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z < 0.2) {
					continue;
				}
				if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				flag_accessible_points_A[i] = true;
				cont_accessible_points_A++;
				cont_unaccessible_points_A--;
			}
		}
		std::cout << "id of orientation:" << ori << std::endl;
		std::cout << "#accessible points of D_A:" << cont_accessible_points_A << std::endl;
		std::cout << "#unaccessible points of D_A:" << cont_unaccessible_points_A << std::endl << std::endl;
	}

//*******************************************//
/////////////////find area S//////////////////
//*******************************************//
	int cont_number_2 = 0;
	std::vector<Eigen::MatrixXd> vis_red_points_DS, vis_red_points_DA;
	std::vector<Eigen::MatrixXd> vis_green_points;
	std::vector<Eigen::MatrixXd> temp_V_vis_DA, temp_V_vis_DS, temp_V_vis_DI;
	for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
		std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I, temp_normal_S, temp_normal_A;
		temp_V_S.resize(sampling_points_in_D_S.size());
		temp_V_A.resize(sampling_points_in_D_A.size());
		temp_V_I.resize(sampling_points_in_D_I.size());
		temp_normal_S.resize(normal_D_S.size());
		temp_normal_A.resize(normal_D_A.size());
		for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
			temp_V_S[i].resize(3, 1);
			temp_V_S[i](0, 0) = sampling_points_in_D_S[i].x();
			temp_V_S[i](1, 0) = sampling_points_in_D_S[i].y();
			temp_V_S[i](2, 0) = sampling_points_in_D_S[i].z();
		}
		for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
			temp_V_A[i].resize(3, 1);
			temp_V_A[i](0, 0) = sampling_points_in_D_A[i].x();
			temp_V_A[i](1, 0) = sampling_points_in_D_A[i].y();
			temp_V_A[i](2, 0) = sampling_points_in_D_A[i].z();
		}
		for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
			temp_V_I[i].resize(3, 1);
			temp_V_I[i](0, 0) = sampling_points_in_D_I[i].x();
			temp_V_I[i](1, 0) = sampling_points_in_D_I[i].y();
			temp_V_I[i](2, 0) = sampling_points_in_D_I[i].z();
		}
		for (int i = 0; i < normal_D_S.size(); i++) {
			temp_normal_S[i].resize(3, 1);
			temp_normal_S[i](0, 0) = normal_D_S[i].x();
			temp_normal_S[i](1, 0) = normal_D_S[i].y();
			temp_normal_S[i](2, 0) = normal_D_S[i].z();
		}
		for (int i = 0; i < normal_D_A.size(); i++) {
			temp_normal_A[i].resize(3, 1);
			temp_normal_A[i](0, 0) = normal_D_A[i].x();
			temp_normal_A[i](1, 0) = normal_D_A[i].y();
			temp_normal_A[i](2, 0) = normal_D_A[i].z();
		}
		temp_V_vis_DS = temp_V_S;
		temp_V_vis_DA = temp_V_A;
		temp_V_vis_DI = temp_V_I;

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < sampling_points_in_D_S.size(); i++)
			temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
		for (int i = 0; i < sampling_points_in_D_A.size(); i++)
			temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
		for (int i = 0; i < sampling_points_in_D_I.size(); i++)
			temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];
		for (int i = 0; i < normal_D_S.size(); i++)
			temp_normal_S[i] = rotMatrix.inverse() * temp_normal_S[i];
		for (int i = 0; i < normal_D_A.size(); i++)
			temp_normal_A[i] = rotMatrix.inverse() * temp_normal_A[i];

		int cont_unaccessible_points = 0;
		for (int i = 0; i < sampling_points_in_D_S.size(); i++) {
			int index_insert;
			if (flag_accessible_points_S[i] == true) {
				continue;
			}

			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_S[i](0, 0) + Cutter.ball_r * temp_normal_S[i](0, 0);
			center_point.m_y = temp_V_S[i](1, 0) + Cutter.ball_r * temp_normal_S[i](1, 0);
			center_point.m_z = temp_V_S[i](2, 0) + Cutter.ball_r * temp_normal_S[i](2, 0);
			Vector3 center_point_2;
			center_point_2.m_x = center_point_2.m_y = center_point_2.m_z = 0;
			center_point_2.m_x = temp_V_S[i](0, 0) + Cutter.ball_r * temp_normal_S[i](0, 0);
			center_point_2.m_y = temp_V_S[i](1, 0) + Cutter.ball_r * temp_normal_S[i](1, 0);
			center_point_2.m_z = temp_V_S[i](2, 0) + Cutter.ball_r * temp_normal_S[i](2, 0);

			if (ori == 0) {
				std::vector<inaccessible_area> temp_vec_inaccessible_area_DS;
				all_the_inaccessible_area_DS.push_back(temp_vec_inaccessible_area_DS);
				map_index_and_DS_points.insert({ cont_unaccessible_points,i });
				map_index_and_DS_points_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;
			}
			else {
				index_insert = map_index_and_DS_points_inv[i];
			}

			if (ori == 0) {
				vis_red_points_DS.push_back(temp_V_vis_DS[i]);
			}

			for (int ii = 0; ii < sampling_points_in_D_I.size(); ii++) {
				bool jud_collision_2 = false;
				//////////////////////////collision detecion////////////////////////////
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r) {
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height && temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_inaccessible_area_DS[all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
					}
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				if (jud_collision_2 == true) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_inaccessible_area_DS[all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
					}
				}
				////////////////////////////////////////////////////////////////////////
			}
		}

		cont_unaccessible_points = 0;
		for (int i = 0; i < sampling_points_in_D_A.size(); i++) {
			int index_insert;
			if (flag_accessible_points_A[i] == true) {
				continue;
			}

			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_A[i](0, 0);
			center_point.m_y = temp_V_A[i](1, 0);
			center_point.m_z = temp_V_A[i](2, 0);

			if (ori == 0) {
				std::vector<inaccessible_area> temp_vec_inaccessible_area_DA;
				all_the_inaccessible_area_DA.push_back(temp_vec_inaccessible_area_DA);
				map_index_and_DA_points.insert({ cont_unaccessible_points,i });
				map_index_and_DA_points_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;
			}
			else {
				index_insert = map_index_and_DA_points_inv[i];
			}

			if (ori == 0) {
				vis_red_points_DA.push_back(temp_V_vis_DA[i]);
			}

			for (int ii = 0; ii < sampling_points_in_D_I.size(); ii++) {
				bool jud_collision_2 = false;
				//////////////////////////collision detecion////////////////////////////
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z <= 0.2) {
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_inaccessible_area_DA[all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
					}
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				if (jud_collision_2 == true) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_inaccessible_area_DA[all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
					}
				}
				////////////////////////////////////////////////////////////////////////
			}
		}
	}
	//*******************************************//

	std::cout<<"Build scalar field..."<<std::endl;
	////////////get the initial scalar field////////////
	//subtractive
	int cont_covering_points = 0;
	for (int i = 0; i < all_the_inaccessible_area_DS.size(); i++) {
		for (int j = 0; j < all_the_inaccessible_area_DS[i].size(); j++) {
			if (flag_covering_points_S[all_the_inaccessible_area_DS[i][j].id_to_point] == false) {
				map_covering_points_and_DS_points.insert({ cont_covering_points,all_the_inaccessible_area_DS[i][j].id_to_point });
				map_covering_points_and_DS_points_inv.insert({ all_the_inaccessible_area_DS[i][j].id_to_point, cont_covering_points });
				std::vector<inaccessible_area> temp_vec_area_s;
				all_the_covering_points_and_S.push_back(temp_vec_area_s);

				inaccessible_area temp_covering_point(i, all_the_inaccessible_area_DS[i][j].id_ori);
				all_the_covering_points_and_S[all_the_covering_points_and_S.size() - 1].push_back(temp_covering_point);
				flag_covering_points_S[all_the_inaccessible_area_DS[i][j].id_to_point] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = map_covering_points_and_DS_points_inv[all_the_inaccessible_area_DS[i][j].id_to_point];
				inaccessible_area temp_covering_point(i, all_the_inaccessible_area_DS[i][j].id_ori);
				all_the_covering_points_and_S[index_insert].push_back(temp_covering_point);
			}
		}
	}
	//这里标量场只考虑了DS，没考虑DA !!
	/*int max_size = -100000000;
	for (int i = 0; i < all_the_covering_points_and_S.size(); i++) {
		max_size = std::max(max_size, int(all_the_covering_points_and_S[i].size()));
		vis_green_points.push_back(temp_V_vis_DI[map_covering_points_and_DS_points[i]]);
	}
	std::vector<double> color_map;
	color_map.resize(all_the_covering_points_and_S.size());
	for (int i = 0; i < all_the_covering_points_and_S.size(); i++)
		color_map[i] = double(all_the_covering_points_and_S[i].size()) / double(max_size);*/
	

	//additive
	cont_covering_points = 0;
	for (int i = 0; i < all_the_inaccessible_area_DA.size(); i++) {
		for (int j = 0; j < all_the_inaccessible_area_DA[i].size(); j++) {
			if (flag_covering_points_A[all_the_inaccessible_area_DA[i][j].id_to_point] == false) {
				map_covering_points_and_DA_points.insert({ cont_covering_points,all_the_inaccessible_area_DA[i][j].id_to_point });
				map_covering_points_and_DA_points_inv.insert({ all_the_inaccessible_area_DA[i][j].id_to_point, cont_covering_points });
				std::vector<inaccessible_area> temp_vec_area_s;
				all_the_covering_points_and_A.push_back(temp_vec_area_s);

				inaccessible_area temp_covering_point(i, all_the_inaccessible_area_DA[i][j].id_ori);
				all_the_covering_points_and_A[all_the_covering_points_and_A.size() - 1].push_back(temp_covering_point);
				flag_covering_points_A[all_the_inaccessible_area_DA[i][j].id_to_point] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = map_covering_points_and_DA_points_inv[all_the_inaccessible_area_DA[i][j].id_to_point];
				inaccessible_area temp_covering_point(i, all_the_inaccessible_area_DA[i][j].id_ori);
				all_the_covering_points_and_A[index_insert].push_back(temp_covering_point);
			}
		}
	}

	for (int i = 0; i < temp_V_vis_DI.size(); i++)
		vis_green_points.push_back(temp_V_vis_DI[i]);
	int max_size_S = -100000000, max_size_A = -100000000;
	for (int i = 0; i < all_the_covering_points_and_S.size(); i++) {
		max_size_S = std::max(max_size_S, int(all_the_covering_points_and_S[i].size()));
	}
	for (int i = 0; i < all_the_covering_points_and_A.size(); i++) {
		max_size_A = std::max(max_size_A, int(all_the_covering_points_and_A[i].size()));
	}
	if (max_size_S == -100000000)
		max_size_S = 0;
	if (max_size_A == -100000000)
		max_size_A = 0;
	color_map.resize(temp_V_vis_DI.size());
	for (int i = 0; i < all_the_covering_points_and_S.size(); i++) {
		int index_V = map_covering_points_and_DS_points[i];
		color_map[index_V] = double(all_the_covering_points_and_S[i].size()) / (double(max_size_S)+ double(max_size_A));
	}
	for (int i = 0; i < all_the_covering_points_and_A.size(); i++) {
		int index_V = map_covering_points_and_DA_points[i];
		color_map[index_V] += double(all_the_covering_points_and_A[i].size()) / (double(max_size_S) + double(max_size_A));
	}


	/////////////////////////////////////////////////////////
	std::cout << "Scalar field has been build !" << std::endl;
	Visual Vis;
	if (open_vis_red_points == true) {
		Vis.creat_red_ball(rootPath+"output\\" + mesh_target + "\\DS", vis_red_points_DS);
		Vis.creat_red_ball_2(rootPath + "output\\" + mesh_target + "\\DA", vis_red_points_DA);
	}
	record_data[2] = vis_red_points_DS.size();
	record_data[3] = vis_red_points_DA.size();
	if (open_vis_green_points == true) {
		Vis.creat_green_ball(rootPath + "output\\" + mesh_target + "\\", vis_green_points, color_map);

		//绿点值归一化
		double max_color_map = -100000000, min_color_map = 100000000;
		for (int i = 0; i < color_map.size(); i++) {
			max_color_map = std::max(max_color_map, color_map[i]);
			min_color_map = std::min(min_color_map, color_map[i]);
		}
		for (int i = 0; i < color_map.size(); i++) {
			if (max_color_map - min_color_map != 0)
				color_map[i] = (color_map[i] - min_color_map) / (max_color_map - min_color_map);
			else
				color_map[i] = 0;
		}
		for (int i = 0; i < vis_green_points.size(); i++)
			current_green_points.push_back(DataPoint(vis_green_points[i](0, 0), vis_green_points[i](1, 0), vis_green_points[i](2, 0), color_map[i]));
	}
		

	std::cout << "Collision Detection done!!" << endl;
	/*std::cout << "&&&time&&& Collision detection: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Subtractive dependency graph: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Green points generation: " << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;*/
	
}

void ReFab::Update_scalar_field_DS_DA(std::vector<int> index_removed_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration)
{
	vector<int> index_V_in_the_remaining_blocks;
	vector<int> current_Di_points;
	vector<bool> current_flag_accessible_points_S = flag_accessible_points_S;
	vector<bool> current_flag_accessible_points_A = flag_accessible_points_A;
	vector<vector<inaccessible_area>> ori_all_the_inaccessible_area_DS = all_the_inaccessible_area_DS;
	vector<vector<inaccessible_area>> ori_all_the_inaccessible_area_DA = all_the_inaccessible_area_DA;
	vector<vector<int>> ori_num_points_of_ori_in_all_the_area_S, ori_num_points_of_ori_in_all_the_area_A;
	ori_num_points_of_ori_in_all_the_area_S.clear(); ori_num_points_of_ori_in_all_the_area_A.clear();
	current_Di_points.clear(); index_V_in_the_remaining_blocks.clear();
	all_the_covering_points_and_A_2 = all_the_covering_points_and_A;
	all_the_covering_points_and_S_2 = all_the_covering_points_and_S;

	//更新不可达点
	ori_num_points_of_ori_in_all_the_area_S.clear(); ori_num_points_of_ori_in_all_the_area_A.clear();
	ori_num_points_of_ori_in_all_the_area_S.resize(ori_all_the_inaccessible_area_DS.size());
	ori_num_points_of_ori_in_all_the_area_A.resize(ori_all_the_inaccessible_area_DA.size());

	//DS
	if (all_the_covering_points_and_S_2.size() != 0) {
		for (int i = 0; i < ori_all_the_inaccessible_area_DS.size(); i++) {
			ori_num_points_of_ori_in_all_the_area_S[i].resize(sample_points_sphere.size());
			for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_S[i].size(); j++)
				ori_num_points_of_ori_in_all_the_area_S[i][j] = 0;
		}
		for (int i = 0; i < ori_all_the_inaccessible_area_DS.size(); i++) {
			for (int itr = 0; itr < ori_all_the_inaccessible_area_DS[i].size(); itr++) {
				ori_num_points_of_ori_in_all_the_area_S[i][ori_all_the_inaccessible_area_DS[i][itr].id_ori]++;
			}
		}
		for (int i = 0; i < index_removed_Di_points.size(); i++) {
			for (int j = 0; j < all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]].size(); j++) {
				int index = all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]][j].id_to_point;
				int ori = all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]][j].id_ori;
				ori_num_points_of_ori_in_all_the_area_S[index][ori]--;
			}
		}
		for (int i = 0; i < ori_num_points_of_ori_in_all_the_area_S.size(); i++) {
			bool flag_inaccessible_points = true;
			for (int j = 0; j < sample_points_sphere.size(); j++) {
				if (ori_num_points_of_ori_in_all_the_area_S[i][j] <= 0) {
					flag_inaccessible_points = false;
					break;
				}
			}
			if (flag_inaccessible_points == false) {
				current_flag_accessible_points_S[map_index_and_DS_points[i]] = true;
				if (flag_accessible_points_S[map_index_and_DS_points[i]] == false) {
					//update all_the_covering_points_and_A since some inaccessible points are removed.
					for (int j = 0; j < all_the_inaccessible_area_DS[i].size(); j++) {
						int index_insert = map_covering_points_and_DS_points_inv[all_the_inaccessible_area_DS[i][j].id_to_point];
						all_the_covering_points_and_S_2[index_insert].erase(all_the_covering_points_and_S_2[index_insert].begin() + all_the_covering_points_and_S_2[index_insert].size() - 1); //这里删除的不是对应的项，只是为了让size - 1
					}
				}
			}
		}

		for (int i = 0; i < current_flag_accessible_points_S.size(); i++) {
			Eigen::MatrixXd temp_mat;
			if (current_flag_accessible_points_S[i] == false) {
				vis_red_points_DS.push_back(temp_mat);
				vis_red_points_DS[vis_red_points_DS.size() - 1].resize(3, 1);
				vis_red_points_DS[vis_red_points_DS.size() - 1](0, 0) = sampling_points_in_D_S[i].x();
				vis_red_points_DS[vis_red_points_DS.size() - 1](1, 0) = sampling_points_in_D_S[i].y();
				vis_red_points_DS[vis_red_points_DS.size() - 1](2, 0) = sampling_points_in_D_S[i].z();
			}
		}
	}

	//DA
	if (all_the_covering_points_and_A_2.size() != 0) {
		for (int i = 0; i < ori_all_the_inaccessible_area_DA.size(); i++) {
			ori_num_points_of_ori_in_all_the_area_A[i].resize(sample_points_sphere.size());
			for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_A[i].size(); j++)
				ori_num_points_of_ori_in_all_the_area_A[i][j] = 0;
		}
		for (int i = 0; i < ori_all_the_inaccessible_area_DA.size(); i++) {
			for (int itr = 0; itr < ori_all_the_inaccessible_area_DA[i].size(); itr++) {
				ori_num_points_of_ori_in_all_the_area_A[i][ori_all_the_inaccessible_area_DA[i][itr].id_ori]++;
			}
		}
		for (int i = 0; i < index_removed_Di_points.size(); i++) {
			for (int j = 0; j < all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]].size(); j++) {
				int index = all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]][j].id_to_point;
				int ori = all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]][j].id_ori;
				ori_num_points_of_ori_in_all_the_area_A[index][ori]--;
			}
		}
		for (int i = 0; i < ori_num_points_of_ori_in_all_the_area_A.size(); i++) {
			bool flag_inaccessible_points = true;
			for (int j = 0; j < sample_points_sphere.size(); j++) {
				if (ori_num_points_of_ori_in_all_the_area_A[i][j] <= 0) {
					flag_inaccessible_points = false;
					break;
				}
			}
			if (flag_inaccessible_points == false) {
				current_flag_accessible_points_A[map_index_and_DA_points[i]] = true;
				if (flag_accessible_points_A[map_index_and_DA_points[i]] == false) {
					//update all_the_covering_points_and_A since some inaccessible points are removed.
					for (int j = 0; j < all_the_inaccessible_area_DA[i].size(); j++) {
						int index_insert = map_covering_points_and_DA_points_inv[all_the_inaccessible_area_DA[i][j].id_to_point];
						all_the_covering_points_and_A_2[index_insert].erase(all_the_covering_points_and_A_2[index_insert].begin() + all_the_covering_points_and_A_2[index_insert].size() - 1); //这里删除的不是对应的项，只是为了让size - 1
					}
				}
			}
		}

		for (int i = 0; i < current_flag_accessible_points_A.size(); i++) {
			Eigen::MatrixXd temp_mat;
			if (current_flag_accessible_points_A[i] == false) {
				vis_red_points_DA.push_back(temp_mat);
				vis_red_points_DA[vis_red_points_DA.size() - 1].resize(3, 1);
				vis_red_points_DA[vis_red_points_DA.size() - 1](0, 0) = sampling_points_in_D_A[i].x();
				vis_red_points_DA[vis_red_points_DA.size() - 1](1, 0) = sampling_points_in_D_A[i].y();
				vis_red_points_DA[vis_red_points_DA.size() - 1](2, 0) = sampling_points_in_D_A[i].z();
			}
		}
	}


}

void ReFab::Update_scalar_field_DS_DA_with_removed_points(std::vector<int> index_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration)
{
	//碰撞检测
	OrientationSamplePoints();
	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_covering_points_S, flag_covering_points_A;
	std::vector<Point_3> sampling_points_from_Di, current_sampling_points_in_Di;
	std::vector<int> index_Di_points_and_current_Di_points;
	index_Di_points_and_current_Di_points.resize(sampling_points_in_D_I.size());
	for (int i = 0; i < index_Di_points.size(); i++) {
		sampling_points_from_Di.push_back(sampling_points_in_D_I[index_Di_points[i]]);
	}
	std::vector<bool> current_flag_accessible_points_S, current_flag_accessible_points_A;
	current_flag_accessible_points_S.resize(sampling_points_from_Di.size());
	current_flag_accessible_points_A.resize(sampling_points_from_Di.size());

	//update current_sampling_points_in_Di by delete index_Di_points from sampling_points_in_D_I
	for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
		bool flag = false;
		for (int j = 0; j < index_Di_points.size(); j++) {
			if (i == index_Di_points[j]) {
				flag = true;
				break;
			}
		}
		if (flag == false) {
			current_sampling_points_in_Di.push_back(sampling_points_in_D_I[i]);
			index_Di_points_and_current_Di_points[i] = (current_sampling_points_in_Di.size() - 1);
		}
	}
	flag_covering_points_S.resize(current_sampling_points_in_Di.size());
	flag_covering_points_A.resize(current_sampling_points_in_Di.size());
	if (current_sampling_points_in_Di.size() == 0)  //新增的
		return;

	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
		current_flag_accessible_points_S[i] = false;
	}
	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
		current_flag_accessible_points_A[i] = false;
	}
	for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
		flag_covering_points_S[i] = false;
		flag_covering_points_A[i] = false;
	}

	//碰撞检测，计算不可达点（任意方向都与Di上的点碰撞）
	for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
		///////////////////rotate/////////////////////
		std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I;
		temp_V_S.resize(sampling_points_from_Di.size());
		temp_V_A.resize(sampling_points_from_Di.size());
		temp_V_I.resize(current_sampling_points_in_Di.size());
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			temp_V_S[i].resize(3, 1);
			temp_V_S[i](0, 0) = sampling_points_from_Di[i].x();
			temp_V_S[i](1, 0) = sampling_points_from_Di[i].y();
			temp_V_S[i](2, 0) = sampling_points_from_Di[i].z();
		}
		temp_V_A = temp_V_S;
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
			temp_V_I[i].resize(3, 1);
			temp_V_I[i](0, 0) = current_sampling_points_in_Di[i].x();
			temp_V_I[i](1, 0) = current_sampling_points_in_Di[i].y();
			temp_V_I[i](2, 0) = current_sampling_points_in_Di[i].z();
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++)
			temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];
		//////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points_A = 0, cont_unaccessible_points_A = 0, cont_accessible_points_S = 0, cont_unaccessible_points_S = 0;
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			if (current_flag_accessible_points_S[i] == true) {
				cont_accessible_points_S++;
				continue;
			}
			cont_unaccessible_points_S++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_S[i](0, 0);
			center_point.m_y = temp_V_S[i](1, 0);
			center_point.m_z = temp_V_S[i](2, 0);
			Vector3 center_point_2;
			center_point_2 = center_point;
			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if ((temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r) || (temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r)) {
					continue;
				}
				if ((temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) && (temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height)) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				current_flag_accessible_points_S[i] = true;
			}
		}

		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			if (current_flag_accessible_points_A[i] == true) {
				cont_accessible_points_A++;
				continue;
			}
			cont_unaccessible_points_A++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_A[i](0, 0);
			center_point.m_y = temp_V_A[i](1, 0);
			center_point.m_z = temp_V_A[i](2, 0);
			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z < 0.2) {
					continue;
				}
				if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				current_flag_accessible_points_A[i] = true;
			}
		}
	}

	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	int cont_number_2 = 0;
	std::vector<std::vector<inaccessible_area>> current_all_the_inaccessible_area_DS, current_all_the_inaccessible_area_DA;
	std::map<int, int> current_map_index_and_DS_points, current_map_index_and_DA_points;
	std::map<int, int> current_map_index_and_DS_points_inv, current_map_index_and_DA_points_inv;
	std::map<int, int> current_map_covering_points_and_DS_points;
	std::map<int, int> current_map_covering_points_and_DS_points_inv;
	std::map<int, int> current_map_covering_points_and_DA_points;
	std::map<int, int> current_map_covering_points_and_DA_points_inv;
	std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_S, current_all_the_covering_points_and_A;
	std::vector<double> current_color_map;
	current_all_the_inaccessible_area_DS.clear();
	std::vector<Eigen::MatrixXd> vis_green_points;
	std::vector<Eigen::MatrixXd> temp_V_vis_DA, temp_V_vis_DS, temp_V_vis_DI;
	for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
		std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I, temp_normal_S, temp_normal_A;
		temp_V_S.resize(sampling_points_from_Di.size());
		temp_V_A.resize(sampling_points_from_Di.size());
		temp_V_I.resize(current_sampling_points_in_Di.size());
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			temp_V_S[i].resize(3, 1);
			temp_V_S[i](0, 0) = sampling_points_from_Di[i].x();
			temp_V_S[i](1, 0) = sampling_points_from_Di[i].y();
			temp_V_S[i](2, 0) = sampling_points_from_Di[i].z();
		}
		temp_V_A = temp_V_S;
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
			temp_V_I[i].resize(3, 1);
			temp_V_I[i](0, 0) = current_sampling_points_in_Di[i].x();
			temp_V_I[i](1, 0) = current_sampling_points_in_Di[i].y();
			temp_V_I[i](2, 0) = current_sampling_points_in_Di[i].z();
		}
		temp_V_vis_DS = temp_V_S;
		temp_V_vis_DA = temp_V_A;
		temp_V_vis_DI = temp_V_I;

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++)
			temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];

		int cont_unaccessible_points = 0;
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			int index_insert;
			if (current_flag_accessible_points_S[i] == true) {
				continue;
			}
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_S[i](0, 0);
			center_point.m_y = temp_V_S[i](1, 0);
			center_point.m_z = temp_V_S[i](2, 0);
			Vector3 center_point_2;
			center_point_2.m_x = center_point_2.m_y = center_point_2.m_z = 0;
			center_point_2.m_x = temp_V_S[i](0, 0);
			center_point_2.m_y = temp_V_S[i](1, 0);
			center_point_2.m_z = temp_V_S[i](2, 0);

			if (ori == 0) {
				std::vector<inaccessible_area> temp_vec_inaccessible_area_DS;
				current_all_the_inaccessible_area_DS.push_back(temp_vec_inaccessible_area_DS);
				current_map_index_and_DS_points.insert({ cont_unaccessible_points,i });
				current_map_index_and_DS_points_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;
			}
			else {
				index_insert = current_map_index_and_DS_points_inv[i];
			}

			if (ori == 0) {
				vis_red_points_DS.push_back(temp_V_vis_DS[i]);
			}

			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				bool jud_collision_2 = false;
				//////////////////////////collision detecion////////////////////////////
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r) {
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height && temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						current_all_the_inaccessible_area_DS[current_all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						current_all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
					}
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				if (jud_collision_2 == true) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						current_all_the_inaccessible_area_DS[current_all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						current_all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
					}
				}
				////////////////////////////////////////////////////////////////////////
			}
		}

		cont_unaccessible_points = 0;
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			int index_insert;
			if (current_flag_accessible_points_A[i] == true) {
				continue;
			}

			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_A[i](0, 0);
			center_point.m_y = temp_V_A[i](1, 0);
			center_point.m_z = temp_V_A[i](2, 0);

			if (ori == 0) {
				std::vector<inaccessible_area> temp_vec_inaccessible_area_DA;
				current_all_the_inaccessible_area_DA.push_back(temp_vec_inaccessible_area_DA);
				current_map_index_and_DA_points.insert({ cont_unaccessible_points,i });
				current_map_index_and_DA_points_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;
			}
			else {
				index_insert = current_map_index_and_DA_points_inv[i];
			}

			if (ori == 0) {
				vis_red_points_DA.push_back(temp_V_vis_DA[i]);
			}

			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				bool jud_collision_2 = false;
				//////////////////////////collision detecion////////////////////////////
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z <= 0.2) {
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						current_all_the_inaccessible_area_DA[current_all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						current_all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
					}
					continue;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
						jud_collision_2 = true;
					}
					else
						continue;
				}
				if (jud_collision_2 == true) {
					inaccessible_area temp_area_S(ii, ori);
					if (ori == 0) {
						current_all_the_inaccessible_area_DA[current_all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
					}
					else {
						cont_number_2++;
						current_all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
					}
				}
				////////////////////////////////////////////////////////////////////////
			}
		}
	}
	//*******************************************//

	std::cout << "Build scalar field..." << std::endl;
	////////////get the initial scalar field////////////
	//subtractive
	int cont_covering_points = 0;
	for (int i = 0; i < current_all_the_inaccessible_area_DS.size(); i++) {
		for (int j = 0; j < current_all_the_inaccessible_area_DS[i].size(); j++) {
			if (flag_covering_points_S[current_all_the_inaccessible_area_DS[i][j].id_to_point] == false) {
				current_map_covering_points_and_DS_points.insert({ cont_covering_points,current_all_the_inaccessible_area_DS[i][j].id_to_point });
				current_map_covering_points_and_DS_points_inv.insert({ current_all_the_inaccessible_area_DS[i][j].id_to_point, cont_covering_points });
				std::vector<inaccessible_area> temp_vec_area_s;
				current_all_the_covering_points_and_S.push_back(temp_vec_area_s);

				inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DS[i][j].id_ori);
				current_all_the_covering_points_and_S[current_all_the_covering_points_and_S.size() - 1].push_back(temp_covering_point);
				flag_covering_points_S[current_all_the_inaccessible_area_DS[i][j].id_to_point] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = current_map_covering_points_and_DS_points_inv[current_all_the_inaccessible_area_DS[i][j].id_to_point];
				inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DS[i][j].id_ori);
				current_all_the_covering_points_and_S[index_insert].push_back(temp_covering_point);
			}
		}
	}
	//additive
	cont_covering_points = 0;
	for (int i = 0; i < current_all_the_inaccessible_area_DA.size(); i++) {
		for (int j = 0; j < current_all_the_inaccessible_area_DA[i].size(); j++) {
			if (flag_covering_points_A[current_all_the_inaccessible_area_DA[i][j].id_to_point] == false) {
				current_map_covering_points_and_DA_points.insert({ cont_covering_points,current_all_the_inaccessible_area_DA[i][j].id_to_point });
				current_map_covering_points_and_DA_points_inv.insert({ current_all_the_inaccessible_area_DA[i][j].id_to_point, cont_covering_points });
				std::vector<inaccessible_area> temp_vec_area_s;
				current_all_the_covering_points_and_A.push_back(temp_vec_area_s);

				inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DA[i][j].id_ori);
				current_all_the_covering_points_and_A[current_all_the_covering_points_and_A.size() - 1].push_back(temp_covering_point);
				flag_covering_points_A[current_all_the_inaccessible_area_DA[i][j].id_to_point] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = current_map_covering_points_and_DA_points_inv[current_all_the_inaccessible_area_DA[i][j].id_to_point];
				inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DA[i][j].id_ori);
				current_all_the_covering_points_and_A[index_insert].push_back(temp_covering_point);
			}
		}
	}

	//更新目前Di与Ds,DA的碰撞关系
	std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_S_2 = all_the_covering_points_and_S_2;
	if (current_all_the_covering_points_and_S_2.size() != 0)
		for (int i = 0; i < index_Di_points.size(); i++) {
			int index = map_covering_points_and_DS_points_inv[index_Di_points[i]];
			current_all_the_covering_points_and_S_2[index].clear();
		}
	std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_A_2 = all_the_covering_points_and_A_2;
	if (current_all_the_covering_points_and_A_2.size() != 0)
		for (int i = 0; i < index_Di_points.size(); i++) {
			int index = map_covering_points_and_DA_points_inv[index_Di_points[i]];
			current_all_the_covering_points_and_A_2[index].clear();
		}
	//更新场
	for (int i = 0; i < temp_V_vis_DI.size(); i++)
		vis_green_points.push_back(temp_V_vis_DI[i]);
	int max_size_S = -100000000, max_size_A = -100000000;
	for (int i = 0; i < current_all_the_covering_points_and_S.size(); i++) {
		max_size_S = std::max(max_size_S, int(current_all_the_covering_points_and_S[i].size()));
	}
	for (int i = 0; i < current_all_the_covering_points_and_S_2.size(); i++) {
		max_size_S = std::max(max_size_S, int(current_all_the_covering_points_and_S_2[i].size()));
	}
	for (int i = 0; i < current_all_the_covering_points_and_A.size(); i++) {
		max_size_A = std::max(max_size_A, int(current_all_the_covering_points_and_A[i].size()));
	}
	for (int i = 0; i < current_all_the_covering_points_and_A_2.size(); i++) {
		max_size_A = std::max(max_size_A, int(current_all_the_covering_points_and_A_2[i].size()));
	}
	current_color_map.resize(temp_V_vis_DI.size());
	for (int i = 0; i < current_color_map.size(); i++) {
		current_color_map[i] = 0;
	}
	for (int i = 0; i < current_all_the_covering_points_and_S.size(); i++) {
		int index_V = current_map_covering_points_and_DS_points[i];
		current_color_map[index_V] = double(current_all_the_covering_points_and_S[i].size()) / (double(max_size_S) + double(max_size_A));
	}
	for (int i = 0; i < current_all_the_covering_points_and_S_2.size(); i++) {
		int index_V = index_Di_points_and_current_Di_points[map_covering_points_and_DS_points[i]];
		current_color_map[index_V] += double(current_all_the_covering_points_and_S_2[i].size()) / (double(max_size_S) + double(max_size_A));
	}
	for (int i = 0; i < current_all_the_covering_points_and_A.size(); i++) {
		int index_V = current_map_covering_points_and_DA_points[i];
		current_color_map[index_V] += double(current_all_the_covering_points_and_A[i].size()) / (double(max_size_S) + double(max_size_A));
	}
	for (int i = 0; i < current_all_the_covering_points_and_A_2.size(); i++) {
		int index_V = index_Di_points_and_current_Di_points[map_covering_points_and_DA_points[i]];
		current_color_map[index_V] += double(current_all_the_covering_points_and_A_2[i].size()) / (double(max_size_S) + double(max_size_A));
	}

	//std::vector<double> current_color_map_2 = current_color_map;
	//for (int i = 0; i < current_color_map.size(); i++) {
	//	current_color_map_2[i] += color_map[index_current_Di_points_and_Di_points[i]];
	//}

	//将切除点加到DS和DA中；可视化
	for (int i = 0; i < current_flag_accessible_points_S.size(); i++) {
		Eigen::MatrixXd temp_mat;
		if (current_flag_accessible_points_S[i] == false) {
			vis_red_points_DS.push_back(temp_mat);
			vis_red_points_DS[vis_red_points_DS.size() - 1].resize(3, 1);
			vis_red_points_DS[vis_red_points_DS.size() - 1](0, 0) = sampling_points_from_Di[i].x();
			vis_red_points_DS[vis_red_points_DS.size() - 1](1, 0) = sampling_points_from_Di[i].y();
			vis_red_points_DS[vis_red_points_DS.size() - 1](2, 0) = sampling_points_from_Di[i].z();
		}
	}
	for (int i = 0; i < current_flag_accessible_points_A.size(); i++) {
		Eigen::MatrixXd temp_mat;
		if (current_flag_accessible_points_A[i] == false) {
			vis_red_points_DA.push_back(temp_mat);
			vis_red_points_DA[vis_red_points_DA.size() - 1].resize(3, 1);
			vis_red_points_DA[vis_red_points_DA.size() - 1](0, 0) = sampling_points_from_Di[i].x();
			vis_red_points_DA[vis_red_points_DA.size() - 1](1, 0) = sampling_points_from_Di[i].y();
			vis_red_points_DA[vis_red_points_DA.size() - 1](2, 0) = sampling_points_from_Di[i].z();
		}
	}
	Visual Vis;
	if (iteration != -1) {
		if (open_vis_red_points == true) {
			Vis.creat_red_ball(rootPath + "temp_vis\\" + mesh_target + "\\New_DS-" + to_string(iteration), vis_red_points_DS);
			Vis.creat_red_ball(rootPath + "temp_vis\\" + mesh_target + "\\New_DA-" + to_string(iteration), vis_red_points_DA);
		}
		/*if (open_vis_green_points == true) {
			Vis.creat_green_ball(rootPath + "temp_vis\\" + mesh_target + "\\New_DI-" + to_string(iteration), vis_green_points, current_color_map);
			current_green_points;
			for (int i = 0; i < vis_green_points.size(); i++)
				current_green_points.push_back(DataPoint(vis_green_points[i](1,0), vis_green_points[i](2, 0), vis_green_points[i](3, 0), current_color_map[i]));
		}*/
	}

	//////暂时测试
	/*ofstream test_interpolation("test_interpolation.txt");
	ofstream test_gradient("gradient.txt");
	std::vector<DataPoint> vertices(vis_green_points.size());
	for (int i = 0; i < vis_green_points.size(); i++) {
		vertices[i].point = Eigen::Vector3d(vis_green_points[i](0, 0), vis_green_points[i](1, 0), vis_green_points[i](2, 0));
		vertices[i].value = current_color_map[i];
	}
	for(double x =-20;x < -15;x += 0.8)
		for(double y = 0;y < 5;y += 0.8)
			for(double z = 65;z < 70;z += 0.8)
			{
				Eigen::Vector3d test_point(x, y, z);
				double interpolatedValue = trilinearInterpolation(test_point, vertices);
				Eigen::Vector3d gradient = gradientField(vertices, test_point, 0.001);
				test_interpolation << test_point.x() << " " << test_point.y() << " " << test_point.z() << " " << interpolatedValue << endl;
				test_gradient << test_point.x() << " " << test_point.y() << " " << test_point.z() << " " << gradient.x() << " " << gradient.y() << " " << gradient.z() << endl;
			}
	cout << "df" << endl;*/

}
//index_removed_Di_points are the index of removed points compared with initial Di.
void ReFab::Update_scalar_field_DS_DA_2(std::vector<int> index_removed_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration)
{
	vector<int> index_V_in_the_remaining_blocks;
	vector<int> current_Di_points;
	vector<bool> current_flag_accessible_points_S = flag_accessible_points_S;
	vector<bool> current_flag_accessible_points_A = flag_accessible_points_A;
	vector<vector<inaccessible_area>> ori_all_the_inaccessible_area_DS = all_the_inaccessible_area_DS;
	vector<vector<inaccessible_area>> ori_all_the_inaccessible_area_DA = all_the_inaccessible_area_DA;
	vector<vector<int>> ori_num_points_of_ori_in_all_the_area_S, ori_num_points_of_ori_in_all_the_area_A;
	ori_num_points_of_ori_in_all_the_area_S.clear(); ori_num_points_of_ori_in_all_the_area_A.clear();
	current_Di_points.clear(); index_V_in_the_remaining_blocks.clear();
	all_the_covering_points_and_A_2 = all_the_covering_points_and_A;
	all_the_covering_points_and_S_2 = all_the_covering_points_and_S;

	//更新不可达点
	ori_num_points_of_ori_in_all_the_area_S.clear(); ori_num_points_of_ori_in_all_the_area_A.clear();
	ori_num_points_of_ori_in_all_the_area_S.resize(ori_all_the_inaccessible_area_DS.size());
	ori_num_points_of_ori_in_all_the_area_A.resize(ori_all_the_inaccessible_area_DA.size());

	//DS
	if (all_the_covering_points_and_S_2.size() != 0) {
		clock_t start_time3 = clock();
		for (int i = 0; i < ori_all_the_inaccessible_area_DS.size(); i++) {
			ori_num_points_of_ori_in_all_the_area_S[i].resize(sample_points_sphere.size());
			for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_S[i].size(); j++)
				ori_num_points_of_ori_in_all_the_area_S[i][j] = 0;
		}
		for (int i = 0; i < ori_all_the_inaccessible_area_DS.size(); i++) {
			for (int itr = 0; itr < ori_all_the_inaccessible_area_DS[i].size(); itr++) {
				ori_num_points_of_ori_in_all_the_area_S[i][ori_all_the_inaccessible_area_DS[i][itr].id_ori]++;
			}
		}
		for (int i = 0; i < index_removed_Di_points.size(); i++) {
			for (int j = 0; j < all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]].size(); j++) {
				int index = all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]][j].id_to_point;
				int ori = all_the_covering_points_and_S[map_covering_points_and_DS_points_inv[index_removed_Di_points[i]]][j].id_ori;
				ori_num_points_of_ori_in_all_the_area_S[index][ori]--;
			}
		}
		clock_t end_time3 = clock();
		//std::cout << "C1 time: " << (double)(end_time3 - start_time3) / CLOCKS_PER_SEC << "s" << std::endl;
		clock_t start_time6 = clock();
		for (int i = 0; i < ori_num_points_of_ori_in_all_the_area_S.size(); i++) {
			bool flag_inaccessible_points = true;
			for (int j = 0; j < sample_points_sphere.size(); j++) {
				if (ori_num_points_of_ori_in_all_the_area_S[i][j] <= 0) {
					flag_inaccessible_points = false;
					break;
				}
			}
			if (flag_inaccessible_points == false) {
				current_flag_accessible_points_S[map_index_and_DS_points[i]] = true;
				//if (flag_accessible_points_S[map_index_and_DS_points[i]] == false) {  //后面似乎用不到，暂时先注释
				//	//update all_the_covering_points_and_A since some inaccessible points are removed.
				//	for (int j = 0; j < all_the_inaccessible_area_DS[i].size(); j++) {
				//		int index_insert = map_covering_points_and_DS_points_inv[all_the_inaccessible_area_DS[i][j].id_to_point];
				//		all_the_covering_points_and_S_2[index_insert].erase(all_the_covering_points_and_S_2[index_insert].begin() + all_the_covering_points_and_S_2[index_insert].size() - 1); //这里删除的不是对应的项，只是为了让size - 1
				//	}
				//}
			}
		}

		for (int i = 0; i < current_flag_accessible_points_S.size(); i++) {
			Eigen::MatrixXd temp_mat;
			if (current_flag_accessible_points_S[i] == false) {
				vis_red_points_DS.push_back(temp_mat);
				vis_red_points_DS[vis_red_points_DS.size() - 1].resize(3, 1);
				vis_red_points_DS[vis_red_points_DS.size() - 1](0, 0) = sampling_points_in_D_S[i].x();
				vis_red_points_DS[vis_red_points_DS.size() - 1](1, 0) = sampling_points_in_D_S[i].y();
				vis_red_points_DS[vis_red_points_DS.size() - 1](2, 0) = sampling_points_in_D_S[i].z();
			}
		}
		clock_t end_time6 = clock();
		//std::cout << "C2 time: " << (double)(end_time6 - start_time6) / CLOCKS_PER_SEC << "s" << std::endl;
	}

	//DA
	if (all_the_covering_points_and_A_2.size() != 0) {
		clock_t start_time3 = clock();
		for (int i = 0; i < ori_all_the_inaccessible_area_DA.size(); i++) {
			ori_num_points_of_ori_in_all_the_area_A[i].resize(sample_points_sphere.size());
			for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_A[i].size(); j++)
				ori_num_points_of_ori_in_all_the_area_A[i][j] = 0;
		}
		for (int i = 0; i < ori_all_the_inaccessible_area_DA.size(); i++) {
			for (int itr = 0; itr < ori_all_the_inaccessible_area_DA[i].size(); itr++) {
				ori_num_points_of_ori_in_all_the_area_A[i][ori_all_the_inaccessible_area_DA[i][itr].id_ori]++;
			}
		}
		clock_t end_time3 = clock();
		clock_t start_time6 = clock();
		for (int i = 0; i < index_removed_Di_points.size(); i++) {
			for (int j = 0; j < all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]].size(); j++) {
				int index = all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]][j].id_to_point;
				int ori = all_the_covering_points_and_A[map_covering_points_and_DA_points_inv[index_removed_Di_points[i]]][j].id_ori;
				ori_num_points_of_ori_in_all_the_area_A[index][ori]--;
			}
		}

		for (int i = 0; i < ori_num_points_of_ori_in_all_the_area_A.size(); i++) {
			bool flag_inaccessible_points = true;
			for (int j = 0; j < sample_points_sphere.size(); j++) {
				if (ori_num_points_of_ori_in_all_the_area_A[i][j] <= 0) {
					flag_inaccessible_points = false;
					break;
				}
			}
			if (flag_inaccessible_points == false) {  
				current_flag_accessible_points_A[map_index_and_DA_points[i]] = true;
				//if (flag_accessible_points_A[map_index_and_DA_points[i]] == false) { //后面似乎用不到，暂时先注释
				//	//update all_the_covering_points_and_A since some inaccessible points are removed.  
				//	for (int j = 0; j < all_the_inaccessible_area_DA[i].size(); j++) {
				//		int index_insert = map_covering_points_and_DA_points_inv[all_the_inaccessible_area_DA[i][j].id_to_point];
				//		all_the_covering_points_and_A_2[index_insert].erase(all_the_covering_points_and_A_2[index_insert].begin() + all_the_covering_points_and_A_2[index_insert].size() - 1); //这里删除的不是对应的项，只是为了让size - 1
				//	}
				//}
			}
		}
		clock_t end_time6 = clock();
		//std::cout << "C4 time: " << (double)(end_time6 - start_time6) / CLOCKS_PER_SEC << "s" << std::endl;
		for (int i = 0; i < current_flag_accessible_points_A.size(); i++) {
			Eigen::MatrixXd temp_mat;
			if (current_flag_accessible_points_A[i] == false) {
				vis_red_points_DA.push_back(temp_mat);
				vis_red_points_DA[vis_red_points_DA.size() - 1].resize(3, 1);
				vis_red_points_DA[vis_red_points_DA.size() - 1](0, 0) = sampling_points_in_D_A[i].x();
				vis_red_points_DA[vis_red_points_DA.size() - 1](1, 0) = sampling_points_in_D_A[i].y();
				vis_red_points_DA[vis_red_points_DA.size() - 1](2, 0) = sampling_points_in_D_A[i].z();
			}
		}
		
	}
}

void ReFab::Update_scalar_field_DS_DA_with_removed_points_2(std::vector<int> index_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration)
{
	//碰撞检测
	OrientationSamplePoints();
	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_covering_points_S, flag_covering_points_A;
	std::vector<Point_3> sampling_points_from_Di,current_sampling_points_in_Di;
	std::vector<int> index_Di_points_and_current_Di_points;
	index_Di_points_and_current_Di_points.resize(sampling_points_in_D_I.size());
	for (int i = 0; i < index_Di_points.size(); i++) {
		sampling_points_from_Di.push_back(sampling_points_in_D_I[index_Di_points[i]]);
	}
	std::vector<bool> current_flag_accessible_points_S, current_flag_accessible_points_A;
	current_flag_accessible_points_S.resize(sampling_points_from_Di.size());
	current_flag_accessible_points_A.resize(sampling_points_from_Di.size());
	
	//update current_sampling_points_in_Di by delete index_Di_points from sampling_points_in_D_I
	for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
		bool flag = false;
		for (int j = 0; j < index_Di_points.size(); j++) {
			if (i == index_Di_points[j]) {
				flag = true;
				break;
			}
		}
		if (flag == false) {
			current_sampling_points_in_Di.push_back(sampling_points_in_D_I[i]);
			index_Di_points_and_current_Di_points[i] = (current_sampling_points_in_Di.size()-1);
		}
	}
	flag_covering_points_S.resize(current_sampling_points_in_Di.size());
	flag_covering_points_A.resize(current_sampling_points_in_Di.size());
	if (current_sampling_points_in_Di.size() == 0)  //新增的
		return;

	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
		current_flag_accessible_points_S[i] = false;
	}
	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
		current_flag_accessible_points_A[i] = false;
	}
	for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
		flag_covering_points_S[i] = false;
		flag_covering_points_A[i] = false;
	}

	//碰撞检测，计算不可达点（任意方向都与Di上的点碰撞）
	for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
		///////////////////rotate/////////////////////
		std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I;
		temp_V_S.resize(sampling_points_from_Di.size());
		temp_V_A.resize(sampling_points_from_Di.size());
		temp_V_I.resize(current_sampling_points_in_Di.size());
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			temp_V_S[i].resize(3, 1);
			temp_V_S[i](0, 0) = sampling_points_from_Di[i].x();
			temp_V_S[i](1, 0) = sampling_points_from_Di[i].y();
			temp_V_S[i](2, 0) = sampling_points_from_Di[i].z();
		}
		temp_V_A = temp_V_S;
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
			temp_V_I[i].resize(3, 1);
			temp_V_I[i](0, 0) = current_sampling_points_in_Di[i].x();
			temp_V_I[i](1, 0) = current_sampling_points_in_Di[i].y();
			temp_V_I[i](2, 0) = current_sampling_points_in_Di[i].z();
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
		for (int i = 0; i < sampling_points_from_Di.size(); i++)
			temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
		for (int i = 0; i < current_sampling_points_in_Di.size(); i++)
			temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];
		//////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points_A = 0, cont_unaccessible_points_A = 0, cont_accessible_points_S = 0, cont_unaccessible_points_S = 0;
		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			if (current_flag_accessible_points_S[i] == true) {
				cont_accessible_points_S++;
				continue;
			}
			cont_unaccessible_points_S++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_S[i](0, 0);
			center_point.m_y = temp_V_S[i](1, 0);
			center_point.m_z = temp_V_S[i](2, 0);
			Vector3 center_point_2;
			center_point_2 = center_point;
			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if ((temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r) || (temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r)) {
					continue;
				}
				if ((temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) && (temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height)) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				current_flag_accessible_points_S[i] = true;
			}
		}

		for (int i = 0; i < sampling_points_from_Di.size(); i++) {
			if (current_flag_accessible_points_A[i] == true) {
				cont_accessible_points_A++;
				continue;
			}
			cont_unaccessible_points_A++;
			bool jud_collision = false;
			Vector3 center_point;
			center_point.m_x = center_point.m_y = center_point.m_z = 0;
			center_point.m_x = temp_V_A[i](0, 0);
			center_point.m_y = temp_V_A[i](1, 0);
			center_point.m_z = temp_V_A[i](2, 0);
			for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
				//用顶点判断//
				Vector3 v_boundary;
				if (temp_V_I[ii](2, 0) - center_point.m_z < 0.2) {
					continue;
				}
				if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
					jud_collision = true;
					break;
				}
				else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				else {
					if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
						jud_collision = true;
						break;
					}
					else
						continue;
				}
				//////////////////

			}
			if (jud_collision == false) {
				current_flag_accessible_points_A[i] = true;
			}
		}
	}

	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	//int cont_number_2 = 0;
	//std::vector<std::vector<inaccessible_area>> current_all_the_inaccessible_area_DS, current_all_the_inaccessible_area_DA;
	//std::map<int, int> current_map_index_and_DS_points, current_map_index_and_DA_points;
	//std::map<int, int> current_map_index_and_DS_points_inv, current_map_index_and_DA_points_inv;
	//std::map<int, int> current_map_covering_points_and_DS_points;
	//std::map<int, int> current_map_covering_points_and_DS_points_inv;
	//std::map<int, int> current_map_covering_points_and_DA_points;
	//std::map<int, int> current_map_covering_points_and_DA_points_inv;
	//std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_S, current_all_the_covering_points_and_A;
	//std::vector<double> current_color_map;
	//current_all_the_inaccessible_area_DS.clear();
	//std::vector<Eigen::MatrixXd> vis_green_points;
	//std::vector<Eigen::MatrixXd> temp_V_vis_DA, temp_V_vis_DS, temp_V_vis_DI;
	//for (int ori = 0; ori < sample_points_sphere.size(); ori++) {
	//	std::vector<Eigen::MatrixXd> temp_V_S, temp_V_A, temp_V_I, temp_normal_S, temp_normal_A;
	//	temp_V_S.resize(sampling_points_from_Di.size());
	//	temp_V_A.resize(sampling_points_from_Di.size());
	//	temp_V_I.resize(current_sampling_points_in_Di.size());
	//	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
	//		temp_V_S[i].resize(3, 1);
	//		temp_V_S[i](0, 0) = sampling_points_from_Di[i].x();
	//		temp_V_S[i](1, 0) = sampling_points_from_Di[i].y();
	//		temp_V_S[i](2, 0) = sampling_points_from_Di[i].z();
	//	}
	//	temp_V_A = temp_V_S;
	//	for (int i = 0; i < current_sampling_points_in_Di.size(); i++) {
	//		temp_V_I[i].resize(3, 1);
	//		temp_V_I[i](0, 0) = current_sampling_points_in_Di[i].x();
	//		temp_V_I[i](1, 0) = current_sampling_points_in_Di[i].y();
	//		temp_V_I[i](2, 0) = current_sampling_points_in_Di[i].z();
	//	}
	//	temp_V_vis_DS = temp_V_S;
	//	temp_V_vis_DA = temp_V_A;
	//	temp_V_vis_DI = temp_V_I;

	//	Eigen::Vector3d vectorBefore(0, 0, 1);
	//	Eigen::Vector3d vectorAfter(sample_points_sphere[ori].x(), sample_points_sphere[ori].y(), sample_points_sphere[ori].z());
	//	rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
	//	for (int i = 0; i < sampling_points_from_Di.size(); i++)
	//		temp_V_S[i] = rotMatrix.inverse() * temp_V_S[i];
	//	for (int i = 0; i < sampling_points_from_Di.size(); i++)
	//		temp_V_A[i] = rotMatrix.inverse() * temp_V_A[i];
	//	for (int i = 0; i < current_sampling_points_in_Di.size(); i++)
	//		temp_V_I[i] = rotMatrix.inverse() * temp_V_I[i];

	//	int cont_unaccessible_points = 0;
	//	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
	//		int index_insert;
	//		if (current_flag_accessible_points_S[i] == true) {
	//			continue;
	//		}
	//		Vector3 center_point;
	//		center_point.m_x = center_point.m_y = center_point.m_z = 0;
	//		center_point.m_x = temp_V_S[i](0, 0);
	//		center_point.m_y = temp_V_S[i](1, 0);
	//		center_point.m_z = temp_V_S[i](2, 0);
	//		Vector3 center_point_2;
	//		center_point_2.m_x = center_point_2.m_y = center_point_2.m_z = 0;
	//		center_point_2.m_x = temp_V_S[i](0, 0);
	//		center_point_2.m_y = temp_V_S[i](1, 0);
	//		center_point_2.m_z = temp_V_S[i](2, 0);

	//		if (ori == 0) {
	//			std::vector<inaccessible_area> temp_vec_inaccessible_area_DS;
	//			current_all_the_inaccessible_area_DS.push_back(temp_vec_inaccessible_area_DS);
	//			current_map_index_and_DS_points.insert({ cont_unaccessible_points,i });
	//			current_map_index_and_DS_points_inv.insert({ i, cont_unaccessible_points });
	//			cont_unaccessible_points++;
	//		}
	//		else {
	//			index_insert = current_map_index_and_DS_points_inv[i];
	//		}

	//		if (ori == 0) {
	//			vis_red_points_DS.push_back(temp_V_vis_DS[i]);
	//		}

	//		for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
	//			bool jud_collision_2 = false;
	//			//////////////////////////collision detecion////////////////////////////
	//			Vector3 v_boundary;
	//			if (temp_V_I[ii](2, 0) - center_point.m_z <= Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z <= Cutter.ball_r) {
	//				continue;
	//			}
	//			else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height && temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r + Cutter.carriage_height) {
	//				inaccessible_area temp_area_S(ii, ori);
	//				if (ori == 0) {
	//					current_all_the_inaccessible_area_DS[current_all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
	//				}
	//				else {
	//					cont_number_2++;
	//					current_all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
	//				}
	//				continue;
	//			}
	//			else if (temp_V_I[ii](2, 0) - center_point.m_z > Cutter.cylinder_height + Cutter.ball_r || temp_V_I[ii](2, 0) - center_point_2.m_z > Cutter.cylinder_height + Cutter.ball_r) {
	//				if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.carriage_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.carriage_r) {
	//					jud_collision_2 = true;
	//				}
	//				else
	//					continue;
	//			}
	//			else {
	//				if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Cutter.ball_r && sqrt(pow(temp_V_I[ii](0, 0) - center_point_2.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point_2.m_y, 2)) <= Cutter.ball_r) {
	//					jud_collision_2 = true;
	//				}
	//				else
	//					continue;
	//			}
	//			if (jud_collision_2 == true) {
	//				inaccessible_area temp_area_S(ii, ori);
	//				if (ori == 0) {
	//					current_all_the_inaccessible_area_DS[current_all_the_inaccessible_area_DS.size() - 1].push_back(temp_area_S);
	//				}
	//				else {
	//					cont_number_2++;
	//					current_all_the_inaccessible_area_DS[index_insert].push_back(temp_area_S);
	//				}
	//			}
	//			////////////////////////////////////////////////////////////////////////
	//		}
	//	}

	//	cont_unaccessible_points = 0;
	//	for (int i = 0; i < sampling_points_from_Di.size(); i++) {
	//		int index_insert;
	//		if (current_flag_accessible_points_A[i] == true) {
	//			continue;
	//		}

	//		Vector3 center_point;
	//		center_point.m_x = center_point.m_y = center_point.m_z = 0;
	//		center_point.m_x = temp_V_A[i](0, 0);
	//		center_point.m_y = temp_V_A[i](1, 0);
	//		center_point.m_z = temp_V_A[i](2, 0);

	//		if (ori == 0) {
	//			std::vector<inaccessible_area> temp_vec_inaccessible_area_DA;
	//			current_all_the_inaccessible_area_DA.push_back(temp_vec_inaccessible_area_DA);
	//			current_map_index_and_DA_points.insert({ cont_unaccessible_points,i });
	//			current_map_index_and_DA_points_inv.insert({ i, cont_unaccessible_points });
	//			cont_unaccessible_points++;
	//		}
	//		else {
	//			index_insert = current_map_index_and_DA_points_inv[i];
	//		}

	//		if (ori == 0) {
	//			vis_red_points_DA.push_back(temp_V_vis_DA[i]);
	//		}

	//		for (int ii = 0; ii < current_sampling_points_in_Di.size(); ii++) {
	//			bool jud_collision_2 = false;
	//			//////////////////////////collision detecion////////////////////////////
	//			Vector3 v_boundary;
	//			if (temp_V_I[ii](2, 0) - center_point.m_z <= 0.2) {
	//				continue;
	//			}
	//			else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle__H_total + 0.2) {
	//				inaccessible_area temp_area_S(ii, ori);
	//				if (ori == 0) {
	//					current_all_the_inaccessible_area_DA[current_all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
	//				}
	//				else {
	//					cont_number_2++;
	//					current_all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
	//				}
	//				continue;
	//			}
	//			else if (temp_V_I[ii](2, 0) - center_point.m_z > Nozzle.nozzle_H_half + 0.2) {
	//				if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.upper_surface_r) {
	//					jud_collision_2 = true;
	//				}
	//				else
	//					continue;
	//			}
	//			else {
	//				if (sqrt(pow(temp_V_I[ii](0, 0) - center_point.m_x, 2) + pow(temp_V_I[ii](1, 0) - center_point.m_y, 2)) <= Nozzle.lowwer_surface_r) {
	//					jud_collision_2 = true;
	//				}
	//				else
	//					continue;
	//			}
	//			if (jud_collision_2 == true) {
	//				inaccessible_area temp_area_S(ii, ori);
	//				if (ori == 0) {
	//					current_all_the_inaccessible_area_DA[current_all_the_inaccessible_area_DA.size() - 1].push_back(temp_area_S);
	//				}
	//				else {
	//					cont_number_2++;
	//					current_all_the_inaccessible_area_DA[index_insert].push_back(temp_area_S);
	//				}
	//			}
	//			////////////////////////////////////////////////////////////////////////
	//		}
	//	}
	//}
	//*******************************************//

	 //后面似乎用不到，暂时先注释
	//std::cout << "Build scalar field..." << std::endl;  
	//////////////get the initial scalar field////////////
	////subtractive
	//int cont_covering_points = 0;
	//for (int i = 0; i < current_all_the_inaccessible_area_DS.size(); i++) {
	//	for (int j = 0; j < current_all_the_inaccessible_area_DS[i].size(); j++) {
	//		if (flag_covering_points_S[current_all_the_inaccessible_area_DS[i][j].id_to_point] == false) {
	//			current_map_covering_points_and_DS_points.insert({ cont_covering_points,current_all_the_inaccessible_area_DS[i][j].id_to_point });
	//			current_map_covering_points_and_DS_points_inv.insert({ current_all_the_inaccessible_area_DS[i][j].id_to_point, cont_covering_points });
	//			std::vector<inaccessible_area> temp_vec_area_s;
	//			current_all_the_covering_points_and_S.push_back(temp_vec_area_s);

	//			inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DS[i][j].id_ori);
	//			current_all_the_covering_points_and_S[current_all_the_covering_points_and_S.size() - 1].push_back(temp_covering_point);
	//			flag_covering_points_S[current_all_the_inaccessible_area_DS[i][j].id_to_point] = true;
	//			cont_covering_points++;
	//		}
	//		else {
	//			int index_insert = current_map_covering_points_and_DS_points_inv[current_all_the_inaccessible_area_DS[i][j].id_to_point];
	//			inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DS[i][j].id_ori);
	//			current_all_the_covering_points_and_S[index_insert].push_back(temp_covering_point);
	//		}
	//	}
	//}
	////additive
	//cont_covering_points = 0;
	//for (int i = 0; i < current_all_the_inaccessible_area_DA.size(); i++) {
	//	for (int j = 0; j < current_all_the_inaccessible_area_DA[i].size(); j++) {
	//		if (flag_covering_points_A[current_all_the_inaccessible_area_DA[i][j].id_to_point] == false) {
	//			current_map_covering_points_and_DA_points.insert({ cont_covering_points,current_all_the_inaccessible_area_DA[i][j].id_to_point });
	//			current_map_covering_points_and_DA_points_inv.insert({ current_all_the_inaccessible_area_DA[i][j].id_to_point, cont_covering_points });
	//			std::vector<inaccessible_area> temp_vec_area_s;
	//			current_all_the_covering_points_and_A.push_back(temp_vec_area_s);

	//			inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DA[i][j].id_ori);
	//			current_all_the_covering_points_and_A[current_all_the_covering_points_and_A.size() - 1].push_back(temp_covering_point);
	//			flag_covering_points_A[current_all_the_inaccessible_area_DA[i][j].id_to_point] = true;
	//			cont_covering_points++;
	//		}
	//		else {
	//			int index_insert = current_map_covering_points_and_DA_points_inv[current_all_the_inaccessible_area_DA[i][j].id_to_point];
	//			inaccessible_area temp_covering_point(i, current_all_the_inaccessible_area_DA[i][j].id_ori);
	//			current_all_the_covering_points_and_A[index_insert].push_back(temp_covering_point);
	//		}
	//	}
	//}

	//更新目前Di与Ds,DA的碰撞关系     //后面似乎用不到，暂时先注释
	//std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_S_2 = all_the_covering_points_and_S_2;
	//if (current_all_the_covering_points_and_S_2.size() != 0)
	//	for (int i = 0; i < index_Di_points.size(); i++) {
	//		int index = map_covering_points_and_DS_points_inv[index_Di_points[i]];
	//		current_all_the_covering_points_and_S_2[index].clear();
	//	}
	//std::vector<std::vector<inaccessible_area>> current_all_the_covering_points_and_A_2 = all_the_covering_points_and_A_2;
	//if (current_all_the_covering_points_and_A_2.size() != 0)
	//	for (int i = 0; i < index_Di_points.size(); i++) {
	//		int index = map_covering_points_and_DA_points_inv[index_Di_points[i]];
	//		current_all_the_covering_points_and_A_2[index].clear();
	//	}
	////更新场
	//for (int i = 0; i < temp_V_vis_DI.size(); i++)
	//	vis_green_points.push_back(temp_V_vis_DI[i]);
	//int max_size_S = -100000000, max_size_A = -100000000;
	//for (int i = 0; i < current_all_the_covering_points_and_S.size(); i++) {
	//	max_size_S = std::max(max_size_S, int(current_all_the_covering_points_and_S[i].size()));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_S_2.size(); i++) {
	//	max_size_S = std::max(max_size_S, int(current_all_the_covering_points_and_S_2[i].size()));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_A.size(); i++) {
	//	max_size_A = std::max(max_size_A, int(current_all_the_covering_points_and_A[i].size()));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_A_2.size(); i++) {
	//	max_size_A = std::max(max_size_A, int(current_all_the_covering_points_and_A_2[i].size()));
	//}
	//current_color_map.resize(temp_V_vis_DI.size());
	//for (int i = 0; i < current_color_map.size(); i++) {
	//	current_color_map[i] = 0;
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_S.size(); i++) {
	//	int index_V = current_map_covering_points_and_DS_points[i];
	//	current_color_map[index_V] = double(current_all_the_covering_points_and_S[i].size()) / (double(max_size_S) + double(max_size_A));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_S_2.size(); i++) {
	//	int index_V = index_Di_points_and_current_Di_points[map_covering_points_and_DS_points[i]]; 
	//	current_color_map[index_V] += double(current_all_the_covering_points_and_S_2[i].size()) / (double(max_size_S) + double(max_size_A));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_A.size(); i++) {
	//	int index_V = current_map_covering_points_and_DA_points[i];
	//	current_color_map[index_V] += double(current_all_the_covering_points_and_A[i].size()) / (double(max_size_S) + double(max_size_A));
	//}
	//for (int i = 0; i < current_all_the_covering_points_and_A_2.size(); i++) {
	//	int index_V = index_Di_points_and_current_Di_points[map_covering_points_and_DA_points[i]];
	//	current_color_map[index_V] += double(current_all_the_covering_points_and_A_2[i].size()) / (double(max_size_S) + double(max_size_A));
	//}
	//
	//将切除点加到DS和DA中；可视化
	for (int i = 0; i < current_flag_accessible_points_S.size(); i++) {
		Eigen::MatrixXd temp_mat;
		if (current_flag_accessible_points_S[i] == false) {
			vis_red_points_DS.push_back(temp_mat);
			vis_red_points_DS[vis_red_points_DS.size() - 1].resize(3, 1);
			vis_red_points_DS[vis_red_points_DS.size() - 1](0, 0) = sampling_points_from_Di[i].x();
			vis_red_points_DS[vis_red_points_DS.size() - 1](1, 0) = sampling_points_from_Di[i].y();
			vis_red_points_DS[vis_red_points_DS.size() - 1](2, 0) = sampling_points_from_Di[i].z();
		}
	}
	for (int i = 0; i < current_flag_accessible_points_A.size(); i++) {
		Eigen::MatrixXd temp_mat;
		if (current_flag_accessible_points_A[i] == false) {
			vis_red_points_DA.push_back(temp_mat);
			vis_red_points_DA[vis_red_points_DA.size() - 1].resize(3, 1);
			vis_red_points_DA[vis_red_points_DA.size() - 1](0, 0) = sampling_points_from_Di[i].x();
			vis_red_points_DA[vis_red_points_DA.size() - 1](1, 0) = sampling_points_from_Di[i].y();
			vis_red_points_DA[vis_red_points_DA.size() - 1](2, 0) = sampling_points_from_Di[i].z();
		}
	}
	Visual Vis;
	if (iteration != -1) {
		if (open_vis_red_points == true) {
			Vis.creat_red_ball(rootPath + "MCTS_temp\\"+mesh_target + "\\New_DS-" + to_string(iteration), vis_red_points_DS);
			Vis.creat_red_ball(rootPath + "MCTS_temp\\"+mesh_target + "\\New_DA-" + to_string(iteration), vis_red_points_DA);
		}
		/*if (open_vis_green_points == true) {
			Vis.creat_green_ball(rootPath + "temp_vis\\" + mesh_target + "\\New_DI-" + to_string(iteration), vis_green_points, current_color_map);
			current_green_points;
			for (int i = 0; i < vis_green_points.size(); i++)
				current_green_points.push_back(DataPoint(vis_green_points[i](1,0), vis_green_points[i](2, 0), vis_green_points[i](3, 0), current_color_map[i]));
		}*/
	}

}


void ReFab::Greedy_one_by_one_delete_points_of_Di()
{
	std::vector<int> index_removed_Di_points;
	std::vector<Point_3> temp_sampling_points_in_D_I = sampling_points_in_D_I;
	std::vector<double> temp_color_map = color_map;
	//sort color_map
	std::vector<int> index_color_map;
	for (int i = 0; i < temp_color_map.size(); i++) {
		index_color_map.push_back(i);
	}
	for (int i = 0; i < temp_color_map.size(); i++) {
		for (int j = i + 1; j < temp_color_map.size(); j++) {
			if (temp_color_map[i] < temp_color_map[j]) {
				double temp_double = temp_color_map[i];
				temp_color_map[i] = temp_color_map[j];
				temp_color_map[j] = temp_double;
				int temp_int = index_color_map[i];
				index_color_map[i] = index_color_map[j];
				index_color_map[j] = temp_int;
			}
		}
	}
	int itr = 1;
	for (int i = 0; i < temp_color_map.size(); i++) {
		index_removed_Di_points.push_back(index_color_map[i]);
		if (i % 100 == 0 && i != 0) {
			std::vector<Eigen::MatrixXd> vis_red_points_DS;
			std::vector<Eigen::MatrixXd> vis_red_points_DA;
			vis_red_points_DS.clear(); vis_red_points_DA.clear();
			
			clock_t start_time, end_time;
			start_time = clock();
			Update_scalar_field_DS_DA(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
			Update_scalar_field_DS_DA_with_removed_points(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
			end_time = clock();
			std::cout << "Update_scalar_field_DS_DA_with_removed_points time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

			itr++;
		}
	}

}

double ReFab::Cut_mesh_by_one_plane(Plane current_plane, double sum_value_of_green)
{
	Slicer_2 slicer;
	slicer.clear();
	slicer.load("output\\"+ mesh_target +"\\D_I.obj");
	slicer.normal[0] = current_plane.normal.x();
	slicer.normal[1] = current_plane.normal.y();
	slicer.normal[2] = current_plane.normal.z();
	slicer.origin[0] = current_plane.origin.x();
	slicer.origin[1] = current_plane.origin.y();
	slicer.origin[2] = current_plane.origin.z();
	slicer.cut();

	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(slicer.normal[0], slicer.normal[1], slicer.normal[2]);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	Slicer_2 temp_slicer = slicer;
	for (int j = 0; j < temp_slicer.positions.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = temp_slicer.positions[j][0];
		temp_V(1, 0) = temp_slicer.positions[j][1];
		temp_V(2, 0) = temp_slicer.positions[j][2];
		temp_V = rotMatrix.inverse() * temp_V;
		temp_slicer.positions[j][0] = temp_V(0, 0);
		temp_slicer.positions[j][1] = temp_V(1, 0);
		temp_slicer.positions[j][2] = temp_V(2, 0);
	}
	Eigen::MatrixXd temp_origin;
	temp_origin.resize(3, 1);
	temp_origin(0, 0) = slicer.origin[0];
	temp_origin(1, 0) = slicer.origin[1];
	temp_origin(2, 0) = slicer.origin[2];
	temp_origin = rotMatrix.inverse() * temp_origin;
	Slicer_2 removed_triangles_from_slicer;
	removed_triangles_from_slicer.clear();
	removed_triangles_from_slicer.positions = slicer.positions;
	for (int j = 0; j < slicer.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if (temp_slicer.positions[slicer.triangles[j][k]][2] - temp_origin(2, 0) > 0.0000001) {
				removed_triangles_from_slicer.triangles.push_back(slicer.triangles[j]);
				slicer.triangles.erase(slicer.triangles.begin() + j);
				j--;
				break;
			}
		}
	}
	//将当前连通分量保留的mesh加入到原始切平面中
	Add_triangles_of_cutting_planes(slicer);
	slicer.save("temp_vis\\"+ mesh_target +"\\result_of_candidate_cutting_planes.obj");
	Mesh temp_mesh;
	CGAL::IO::read_OBJ("temp_vis\\"+ mesh_target +"\\result_of_candidate_cutting_planes.obj", temp_mesh);
	double cost_volume =  1- CGAL::Polygon_mesh_processing::volume(temp_mesh) / volume_of_Di;

	////用离散点近似绿厂的积分
	//double current_value_of_green = 0;
	//int cont_green_points = 0;
	//CGAL::Side_of_triangle_mesh<Mesh, K> inside(temp_mesh);
	//for (int i = 0; i < current_green_points.size(); i++) {
	//	Point_3 temp_green_point (current_green_points[i].point.x(), current_green_points[i].point.y(), current_green_points[i].point.z());
	//	CGAL::Bounded_side res = inside(temp_green_point);
	//	if (res != CGAL::ON_BOUNDED_SIDE) {    //切除部分绿点平均值越大越好
	//		current_value_of_green += current_green_points[i].value;
	//		cont_green_points++;
	//	}
	//}

	double cost =  cost_volume;
	return cost;
}

vector<Plane> ReFab::Sort_cutting_planes(vector<Plane> Selected_Planes, std::vector<std::vector<Point_3>>& additive_interface_points, std::vector<std::vector<int>>& coverd_additive_components_by_selected_planes)
{
	Element* element = new Element [Selected_Planes.size()];
	vector<vector<Point_3>> all_the_points_of_cutting_planes (Selected_Planes.size());
	for (int i = 0; i < Selected_Planes.size(); i++) {
		element[i].name = to_string(i);
		Slicer_2 slicer;
		slicer.clear();
		slicer.load("output\\" + mesh_target + "\\D_I.obj");
		slicer.normal[0] = Selected_Planes[i].normal.x();
		slicer.normal[1] = Selected_Planes[i].normal.y();
		slicer.normal[2] = Selected_Planes[i].normal.z();
		slicer.origin[0] = Selected_Planes[i].origin.x();
		slicer.origin[1] = Selected_Planes[i].origin.y();
		slicer.origin[2] = Selected_Planes[i].origin.z();
		slicer.cut();
		std::map<std::pair<int, int>, int>::iterator it;
		for (it = slicer.intersections.begin(); it != slicer.intersections.end(); it++)
			all_the_points_of_cutting_planes[i].push_back(Point_3(slicer.positions[it->second][0], slicer.positions[it->second][1], slicer.positions[it->second][2]));
	}
	/*ofstream off("sdfsdf.obj");
	for (int i = 0; i < all_the_points_of_cutting_planes.size(); i++)
		for (int j = 0; j < all_the_points_of_cutting_planes[i].size(); j++)
			off << "v " << all_the_points_of_cutting_planes[i][j].x() << " " << all_the_points_of_cutting_planes[i][j].y() << " " << all_the_points_of_cutting_planes[i][j].z() << endl;*/

	for (int i = 0; i < Selected_Planes.size(); i++) {
		double a = Selected_Planes[i].normal.x();
		double b = Selected_Planes[i].normal.y();
		double c = Selected_Planes[i].normal.z();
		double d = -a * Selected_Planes[i].origin.x() - b * Selected_Planes[i].origin.y() - c * Selected_Planes[i].origin.z();
		for (int j = 0; j < Selected_Planes.size(); j++) {
			if(i == j)
				continue;
			bool jud_contain = true;
			for (int k = 0; k < all_the_points_of_cutting_planes[j].size(); k++) {
				double distance = (a * all_the_points_of_cutting_planes[j][k].x() + b * all_the_points_of_cutting_planes[j][k].y() + c * all_the_points_of_cutting_planes[j][k].z() + d) / std::sqrt(a * a + b * b + c * c);
				if (distance < 0) {
					jud_contain = false;
					break;
				}
			}
			if (jud_contain == true) {
				element[j].containedBy.insert(&element[i]);
				cout << i << " contains " << j << endl;
			}
		}
	}

	vector<Element*> elements;
	for (int i = 0; i < Selected_Planes.size(); i++)
		elements.push_back(&element[i]);
	Element temp;
	vector<Element*> sorted = temp.topologicalSort(elements);
	for (Element* elem : sorted) {
		cout << elem->name << " "<<endl;
	}
	vector<Plane> New_Selected_Planes(Selected_Planes.size());
	std::vector<std::vector<int>> new_coverd_additive_components_by_selected_planes = coverd_additive_components_by_selected_planes;
	//std::vector<std::vector<Point_3>> new_additive_interface_points = additive_interface_points;
	for (int i = 0; i < sorted.size(); i++) {
		New_Selected_Planes[i] = Selected_Planes[stoi(sorted[i]->name)];
		//new_additive_interface_points[i] = additive_interface_points[stoi(sorted[i]->name)];
		new_coverd_additive_components_by_selected_planes[i] = coverd_additive_components_by_selected_planes[stoi(sorted[i]->name)];
	}
	//additive_interface_points = new_additive_interface_points;
	coverd_additive_components_by_selected_planes = new_coverd_additive_components_by_selected_planes;
	return New_Selected_Planes;
}

void ReFab::solveWeightedSetCover(glp_prob* lp, int n, int m, std::vector<int> weights, std::vector<std::set<int>> subsets)
{
	glp_set_prob_name(lp, "weighted set cover");
	glp_set_obj_dir(lp, GLP_MIN); // 目标是最小化

	// 添加列（决策变量）
	glp_add_cols(lp, m);
	for (int i = 1; i <= m; i++) {
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0); // 0-1 变量
		glp_set_obj_coef(lp, i, weights[i - 1]); // 设置目标函数系数
		glp_set_col_kind(lp, i, GLP_BV); // 设置为二元变量
	}

	// 添加行（约束条件）
	glp_add_rows(lp, n);
	for (int k = 1; k <= n; k++) {
		glp_set_row_bnds(lp, k, GLP_LO, 1.0, 0.0); // 每个元素至少被一个子集覆盖
	}

	// 创建非零元素的数组
	std::vector<int> ia(1 + n * m);
	std::vector<int> ja(1 + n * m);
	std::vector<double> ar(1 + n * m);

	int idx = 1;
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < m; j++) {
			if (subsets[j].count(k + 1) > 0) {
				ia[idx] = k + 1; // 行索引
				ja[idx] = j + 1; // 列索引
				ar[idx] = 1.0; // 系数
				idx++;
			}
		}
	}

	// 加载矩阵
	glp_load_matrix(lp, idx - 1, ia.data(), ja.data(), ar.data());
}

void ReFab::Initialize_Cutting_Planes(string additive_interface)
{
	double arfa = 0.9; 
	float step = 30.0; //offset步长  30

	//Calculate_mean_curvature("models\\kitten.off", true);

	//typedef CGAL::Simple_cartesian<double> Kernel;
	//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
	////load .off file, and then translate to .xyz file
	//std::ifstream input(additive_interface);
	//// 创建Polyhedron对象并读取.off文件
	//Polyhedron polyhedron;
	//input >> polyhedron;
	
	///////////////calculate connected components////////////////////
	additive_interface_points.clear();
	Mesh mesh;
	if (!PMP::IO::read_polygon_mesh(additive_interface, mesh))
	{
		std::cerr << "Invalid input." << std::endl;
	}
	vector<bool> flag_face_been_search(mesh.number_of_faces());
	for (int i = 0; i < mesh.number_of_faces(); i++)
		flag_face_been_search[i] = false;
	typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
	const double bound = std::cos(0.75 * CGAL_PI);
	std::vector<face_descriptor> cc;
	for (face_descriptor fd : faces(mesh))
	{
		cc.clear();
		if (flag_face_been_search[fd.idx()] == true)
			continue;
		PMP::connected_component(fd,
			mesh,
			std::back_inserter(cc));
		std::vector<Point_3> temp_points;
		additive_interface_points.push_back(temp_points);
		for (int j = 0; j < cc.size(); j++) {
			flag_face_been_search[cc[j].idx()] = true;  
			for (Mesh::Vertex_index vd : mesh.vertices_around_face(mesh.halfedge(cc[j]))) {
				Point_3 p = mesh.point(vd);
				additive_interface_points[additive_interface_points.size() - 1].push_back(p);
			}
		}
	}
	/*ofstream ofile("output/aaa.obj");
	for (int i = 0; i < additive_interface_points.size(); i++)
		for (int j = 0; j < additive_interface_points[i].size(); j++)
			ofile << "v " << additive_interface_points[i][j].x() << " " << additive_interface_points[i][j].y() << " " << additive_interface_points[i][j].z() << endl;*/

	std::cout << "number of connected components:" << additive_interface_points.size() << endl;
	/////////////////////////////////////////////////////////////////

	/*std::vector<Point_3> points;
	for (Polyhedron::Vertex_iterator iter = polyhedron.vertices_begin(); iter != polyhedron.vertices_end(); ++iter) {
		Point_3 point(iter->point().x(), iter->point().y(), iter->point().z());
		points.push_back(point);
	}*/


	/////////////////////////计算初始候选切平面/////////////////////////
	//用基底点的平均值作为基底中心
	Vec3 temp_point(0, 0, 0);
	for (int i = 0; i < vertices_base.size(); i++) {
		temp_point.m_x += vertices_base[i].m_x;
		temp_point.m_y += vertices_base[i].m_y;
		temp_point.m_z += vertices_base[i].m_z;
	}
	temp_point /= vertices_base.size();
	Point_3 origin_point_in_base(temp_point.m_x, temp_point.m_y, temp_point.m_z); 

	std::vector<Point_3> sample_normals = OrientationSamplePoints_normal();
	std::vector<std::vector<Plane>> Planes;
	std::vector<std::vector<std::vector<double>>> Cost_planes;  //对应的连通分量[采样方向]
	for (int p = 0; p < additive_interface_points.size(); p++) {
		std::vector<Plane> temp_planes;
		Planes.push_back(temp_planes);
		Cost_planes.push_back(std::vector<std::vector<double>>());
		for (int i = 0; i < sample_normals.size(); i++) {
			bool jud_delete = false;
			double a = sample_normals[i].x();
			double b = sample_normals[i].y();
			double c = sample_normals[i].z();
			double d = -a * origin_point_in_base.x() - b * origin_point_in_base.y() - c * origin_point_in_base.z();
			double min_distance = 100000000;
			int index_closest_point = -1;
			for (int j = 0; j < additive_interface_points[p].size(); j++) {
				double distance = (a * additive_interface_points[p][j].x() + b * additive_interface_points[p][j].y() + c * additive_interface_points[p][j].z() + d) / std::sqrt(a * a + b * b + c * c);
				if (distance < -0.01) { 
					jud_delete = true;
					break;
				}
				if (min_distance > distance) {
					min_distance = distance;
					index_closest_point = j;
				}
			}

			//若存在点在平面下方，则删除该平面
			if (jud_delete == true)
				continue;
			Plane temp_plane;
			Planes[p].push_back(temp_plane);
			Cost_planes[p].push_back(std::vector<double>());
			Planes[p][Planes[p].size() - 1].normal = sample_normals[i];
			Planes[p][Planes[p].size() - 1].origin = additive_interface_points[p][index_closest_point];
			//若基底不完全在平面下方，则删除该平面
			bool jud_continue = false;
			d = -a * additive_interface_points[p][index_closest_point].x() - b * additive_interface_points[p][index_closest_point].y() - c * additive_interface_points[p][index_closest_point].z();
			for (int j = 0; j < vertices_base.size(); j++) {
				double distance = (a * vertices_base[j].m_x + b * vertices_base[j].m_y + c * vertices_base[j].m_z + d) / std::sqrt(a * a + b * b + c * c);
				if (distance >= 0) {
					Planes[p].erase(Planes[p].end()-1);
					Cost_planes[p].erase(Cost_planes[p].end() - 1);
					jud_continue = true;
					break;
				}
			}
			if (jud_continue == true)
				continue;
			for (int t = 0; t < additive_interface_points.size(); t++) {
				bool jud_contain = true;
				for (int j = 0; j < additive_interface_points[t].size(); j++) {
					double distance = (a * additive_interface_points[t][j].x() + b * additive_interface_points[t][j].y() + c * additive_interface_points[t][j].z() + d) / std::sqrt(a * a + b * b + c * c);
					if (distance < 0) {
						jud_contain = false;
						break;
						
					}
				}
				if(jud_contain == false)
					Cost_planes[p][Cost_planes[p].size() - 1].push_back(-1);
				else
					Cost_planes[p][Cost_planes[p].size() - 1].push_back(0);
			}
		}
	} 
	///////////////////////////////////////////////////////////////////////////

	/////////////////////////打分函数挑选初始切平面/////////////////////////
	//用切割后剩余体积作为打分函数分数
	std::cout << "Selecting cutting planes from candidate planes..." << endl;
	Eigen::MatrixXd V_f;
	Eigen::MatrixXi F_f;
	igl::readOFF("output\\"+ mesh_target +"\\D_I.off", V_f, F_f);
	igl::writeOBJ("output\\"+ mesh_target +"\\D_I.obj", V_f, F_f);
	//Selected_Planes.resize(Planes.size());
	Mesh temp_mesh;
	CGAL::IO::read_OFF("output\\" + mesh_target + "\\D_I.off", temp_mesh);
	volume_of_Di = CGAL::Polygon_mesh_processing::volume(temp_mesh);
	double sum_value_of_green = 0;
	for (int i = 0; i < current_green_points.size(); i++)
		sum_value_of_green += current_green_points[i].value;

	for (int i = 0; i < Planes.size(); i++) {
		vector<double> cost(Planes[i].size());
		for (int j = 0; j < Planes[i].size(); j++) {
			//cost[j] = 0;
			//double a = Planes[i][j].normal.x();
			//double b = Planes[i][j].normal.y();
			//double c = Planes[i][j].normal.z();
			//double d = -a * Planes[i][j].origin.x() - b * Planes[i][j].origin.y() - c * Planes[i][j].origin.z();
			//for (int k = 0; k < additive_interface_points[i].size(); k++) {
			//	double distance = a * additive_interface_points[i][k].x() + b * additive_interface_points[i][k].y() + c * additive_interface_points[i][k].z() + d / std::sqrt(a * a + b * b + c * c);
			//	/*if (distance < 0)
			//		cout << "a" << endl;*/
			//	cost[j] += distance;
			//}
			//cost[j] = cost[j] / additive_interface_points[i].size();
			//cost[j] = double(1) / cost[j];
			cost[j] = Cut_mesh_by_one_plane(Planes[i][j], sum_value_of_green);
			Cost_planes[i][j][i] = cost[j];
		}
		//find index from max score
		int index_max_score = -1;
		double max_score = -100000000;
		for (int j = 0; j < cost.size(); j++) {
			if (max_score < cost[j]) {
				max_score = cost[j];
				index_max_score = j;
			}
		}
		//Selected_Planes[i] = Planes[i][index_max_score];
	}
	//build the sets
	vector<bool> is_covered(additive_interface_points.size());
	std::vector<std::set<int>> all_candidate_sets;
	std::vector<int> all_candidate_weights;
	//std::vector<std::set<int>> all_candidate_subsets;
	for (int i = 0; i < Cost_planes.size(); i++) {
		is_covered[i] = false;
		for (int j = 0; j < Cost_planes[i].size(); j++) {
			std::set<int> temp_set;
			for (int k = 0; k < Cost_planes[i][j].size(); k++) {
				if (Cost_planes[i][j][k] != -1) {
					temp_set.insert(k+1);
					if (Cost_planes[i][j][k] != 0)
						all_candidate_weights.push_back(Cost_planes[i][j][k]*1000);
				}
			}
			all_candidate_sets.push_back(temp_set);
		}
	}

	//////////////////////////整数线性规划/////////////////////////////
	// 初始化GLPK问题
	glp_prob* lp = glp_create_prob();
	solveWeightedSetCover(lp, Cost_planes.size(), all_candidate_sets.size(), all_candidate_weights, all_candidate_sets);
	// 求解整数线性规划问题
	glp_simplex(lp, nullptr);
	glp_intopt(lp, nullptr);
	// 输出结果
	std::cout << "Optimal solution found:" << std::endl;
	int cont_set = 0;
	for (int i = 0; i < Planes.size(); i++) {
		for (int j = 0; j < Planes[i].size(); j++) {
			if (glp_mip_col_val(lp, cont_set + 1)==1) {
				Selected_Planes.push_back(Planes[i][j]);
			}
			cont_set++;
		}
	}
	vector<int> save_selected_planes;
	for (int i = 1; i <= all_candidate_sets.size(); i++) {
		std::cout << "x[" << i << "] = " << glp_mip_col_val(lp, i) << std::endl;
		if(glp_mip_col_val(lp, i) == 1)
			save_selected_planes.push_back(i-1);
	}
	// 清理
	glp_delete_prob(lp);
	glp_free_env();
	////////////////////////////////////////////////////////////////////

	//Greedy method to solve Weighted Set Cover Problem
	//bool jud_all_covered;
	//int temp_cont = 0;
	//double sum_cost = 0;
	//coverd_additive_components_by_selected_planes.clear();
	//vector<bool> flag_candidate_set_be_used(all_candidate_sets.size());
	//for(int i = 0; i < all_candidate_sets.size(); i++)
	//	flag_candidate_set_be_used[i] = false;
	//while (true) {
	//	jud_all_covered = true;
	//	double max_score = -100000000;
	//	int index_max_score = -1;
	//	for (int i = 0; i < all_candidate_sets.size(); i++) {
	//		int cont_uncovered = 0;
	//		for (int j = 0; j < all_candidate_sets[i].size(); j++)
	//			if (is_covered[all_candidate_sets[i][j]] == false)
	//				cont_uncovered++;
	//		//double current_score = arfa * (1 / all_candidate_weights[i]) + (1 - arfa) * cont_uncovered; 
	//		double current_score = arfa * (1 - all_candidate_weights[i]) + (1 - arfa) * cont_uncovered / is_covered.size();
	//		if(current_score > max_score && cont_uncovered != 0 && flag_candidate_set_be_used[i] == false) {
	//			max_score = current_score;
	//			index_max_score = i;
	//		}
	//	}
	//	coverd_additive_components_by_selected_planes.push_back(std::vector<int>());
	//	for (int j = 0; j < all_candidate_sets[index_max_score].size(); j++) {
	//		if (is_covered[all_candidate_sets[index_max_score][j]] == false) {
	//			is_covered[all_candidate_sets[index_max_score][j]] = true;
	//			coverd_additive_components_by_selected_planes[coverd_additive_components_by_selected_planes.size() - 1].push_back(all_candidate_sets[index_max_score][j]);
	//		}
	//	}
	//	flag_candidate_set_be_used[index_max_score] = true;
	//	//all_candidate_sets.erase(all_candidate_sets.begin() + index_max_score);
	//	for (int j = 0; j < is_covered.size(); j++) {
	//		if (is_covered[j] == false) {
	//			jud_all_covered = false;
	//			break;
	//		}
	//	}
	//	int cont = 0;
	//	bool jud_find = false;
	//	for (int i = 0; i < Planes.size(); i++) {
	//		for (int j = 0; j < Planes[i].size(); j++) {
	//			if (cont == index_max_score) {
	//				jud_find = true;
	//				Selected_Planes.push_back(Planes[i][j]);
	//				cout << i << j << endl;
	//				break;
	//			}
	//			cont++;
	//		}
	//		if(jud_find == true)
	//			break;
	//	}
	//	temp_cont++;
	//	sum_cost += max_score;
	//	if (jud_all_covered == true) {
	//		break;
	//	}
	//}
	//std::cout << temp_cont << endl;


	coverd_additive_components_by_selected_planes.clear();
	for (int i = 0; i < save_selected_planes.size(); i++) {
		coverd_additive_components_by_selected_planes.push_back(std::vector<int>());
		for (int j = 0; j < all_candidate_sets[save_selected_planes[i]].size(); j++) {
			std::set<int>::iterator it = all_candidate_sets[save_selected_planes[i]].begin();
			std::advance(it, j);
			if (is_covered[*it] == false) {
				is_covered[*it] = true;
				coverd_additive_components_by_selected_planes[coverd_additive_components_by_selected_planes.size() - 1].push_back(*it - 1);
			}
		}
	}
	
	//排序，使每次都能切到mesh
	Selected_Planes = Sort_cutting_planes(Selected_Planes, additive_interface_points, coverd_additive_components_by_selected_planes);
	std::cout << "Selecting done." << endl;

	vector<Eigen::Vector3d> lines; lines.resize(Selected_Planes.size());
	std::vector<Eigen::MatrixXd> vis_points; vis_points.resize(Selected_Planes.size());
	for (int i = 0; i < Selected_Planes.size(); i++) {
		lines[i].x() = Selected_Planes[i].normal.x();
		lines[i].y() = Selected_Planes[i].normal.y();
		lines[i].z() = Selected_Planes[i].normal.z();
		vis_points[i].resize(3, 1);
		vis_points[i](0, 0) = Selected_Planes[i].origin.x();
		vis_points[i](1, 0) = Selected_Planes[i].origin.y();
		vis_points[i](2, 0) = Selected_Planes[i].origin.z();
	}
	string file_name = mesh_target + "\\initial_cutting_planes_normal";
	Visual vis;
	vis.generateModelForRendering(lines, file_name, vis_points);
	///////////////////////////////////////////////////////////////////////////


	/////////////////////////切平面与mesh相交/////////////////////////
	flag_useful_cutting_planes = new bool[Selected_Planes.size()];
	Slicer_2 slicer;
	vector<vector<Vec3>> vis_cutting_planes(Selected_Planes.size());
	Eigen::MatrixXd V_ff;
	Eigen::MatrixXi F_ff;
	igl::readOFF("output\\"+ mesh_target +"\\D_I.off", V_ff, F_ff);
	igl::writeOBJ("output\\"+ mesh_target +"\\D_I.obj", V_ff, F_ff);
	flag_use_selected_planes.resize(Selected_Planes.size());
	//计算基底到切平面的最近距离
	vector<double> min_distance_base_to_plane(Selected_Planes.size());
	for (int i = 0; i < Selected_Planes.size(); i++) {
		double a = Selected_Planes[i].normal.x();
		double b = Selected_Planes[i].normal.y();
		double c = Selected_Planes[i].normal.z();
		double d = -a * Selected_Planes[i].origin.x() - b * Selected_Planes[i].origin.y() - c * Selected_Planes[i].origin.z();
		min_distance_base_to_plane[i] = 100000000;
		for (int j = 0; j < vertices_base.size(); j++) {
			double distance = abs((a * vertices_base[j].m_x + b * vertices_base[j].m_y + c * vertices_base[j].m_z + d)) / std::sqrt(a * a + b * b + c * c);
			if (distance < min_distance_base_to_plane[i])
				min_distance_base_to_plane[i] = distance;
		}
	}

	std::vector<Plane> ori_Selected_Planes = Selected_Planes;
	bool jud_is_force_offset = true;
	for (float percent_offset = 0; percent_offset <= 100; percent_offset+= step){ 
		for (int i = 0; i < Selected_Planes.size(); i++)
			flag_useful_cutting_planes[i] = true;
		slicer.clear();
		int index_last_useful_planes;
		if(percent_offset == 0 || jud_is_force_offset)
			slicer.load("output\\"+ mesh_target +"\\D_I.obj");
		else
			slicer.load("temp_vis\\" + mesh_target + "\\(C)new_D_I-" + to_string(index_last_useful_planes) + ".obj");  //改为只读取上一次offset截取的部分
		//slicer.load("output\\"+ mesh_target +"\\D_I.obj");
		vector<int> index_exist_planes;
		index_exist_planes.clear();
		barycenter_of_cutting_planes.clear();
		for (int i = 0; i < Selected_Planes.size(); i++) {
			//更新切平面
			Selected_Planes[i].origin = Point_3(ori_Selected_Planes[i].origin.x() - ori_Selected_Planes[i].normal.x() * min_distance_base_to_plane[i] * percent_offset / 100.0,
				ori_Selected_Planes[i].origin.y() - ori_Selected_Planes[i].normal.y() * min_distance_base_to_plane[i] * percent_offset / 100.0, 
				ori_Selected_Planes[i].origin.z() - ori_Selected_Planes[i].normal.z() * min_distance_base_to_plane[i] * percent_offset / 100.0);
			if (Selected_Planes[i].normal.x() == 0 && Selected_Planes[i].normal.y() == 0 && Selected_Planes[i].normal.z() == 1)
				Selected_Planes[i].origin = Point_3(Selected_Planes[i].origin.x(), Selected_Planes[i].origin.y(), Selected_Planes[i].origin.z()- 1);
			//若注释该部分，则初始化的切平面数量总是等于增材接触面的连通分量数***********************
			//判断该连通红面是否都已经在已有切面之外
			/*bool jud_all_points_outside = true;
			flag_use_selected_planes[i] = true;
			if (i != 0) {
				for (int p = 0; p < index_exist_planes.size(); p++) {
					double a = Selected_Planes[index_exist_planes[p]].normal.x();
					double b = Selected_Planes[index_exist_planes[p]].normal.y();
					double c = Selected_Planes[index_exist_planes[p]].normal.z();
					double d = -a * Selected_Planes[index_exist_planes[p]].origin.x() - b * Selected_Planes[index_exist_planes[p]].origin.y() - c * Selected_Planes[index_exist_planes[p]].origin.z();
					for (int k = 0; k < additive_interface_points[i].size(); k++) {
						double distance = a * additive_interface_points[i][k].x() + b * additive_interface_points[i][k].y() + c * additive_interface_points[i][k].z() + d / std::sqrt(a * a + b * b + c * c);
						if (distance < 0) {
							jud_all_points_outside = false;
							break;
						}
					}
					if(jud_all_points_outside == false)
						break;
				}
				if (jud_all_points_outside == true) {
					flag_use_selected_planes[i] = false;
					continue;
				}
			}
			index_exist_planes.push_back(i);*/

			Slicer_2 removed_triangles_from_slicer;
			vector<vector<int>> cutting_plane_polygons;
			int ori_num_vertex = slicer.positions.size();
			slicer.normal[0] = Selected_Planes[i].normal.x();
			slicer.normal[1] = Selected_Planes[i].normal.y();
			slicer.normal[2] = Selected_Planes[i].normal.z();
			slicer.origin[0] = Selected_Planes[i].origin.x();
			slicer.origin[1] = Selected_Planes[i].origin.y();
			slicer.origin[2] = Selected_Planes[i].origin.z();
			slicer.cut();
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(Selected_Planes[i].normal.x(), Selected_Planes[i].normal.y(), Selected_Planes[i].normal.z());
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
			temp_cutting_point(0, 0) = Selected_Planes[i].origin.x();
			temp_cutting_point(1, 0) = Selected_Planes[i].origin.y();
			temp_cutting_point(2, 0) = Selected_Planes[i].origin.z();
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

			if (cutting_plane_edges.size() == 0) {
				flag_useful_cutting_planes[i] = false;
				continue;
			}
				
			//get a polygon from the cutting_plane_edges,this polygon is made up of pairs in cutting_plane_edges,the adjacent points in a polygon are always in some pair
			vector<bool> flag_edge_been_search(cutting_plane_edges.size());
			for (int j = 0; j < cutting_plane_edges.size(); j++)
				flag_edge_been_search[j] = false;
			Point_3 temp_point = Selected_Planes[i].origin;
			double min_distance = 100000000;
			int index_closest_point = -1;
			for (int j = 0; j < cutting_plane_edges.size(); j++) {
				double distance = std::sqrt(pow(temp_point.x() - slicer.positions[cutting_plane_edges[j].first][0], 2) + pow(temp_point.y() - slicer.positions[cutting_plane_edges[j].first][1], 2) + pow(temp_point.z() - slicer.positions[cutting_plane_edges[j].first][2], 2));
				if (min_distance > distance) {
					min_distance = distance;
					index_closest_point = j;
				}
			}
			/*if (flag_edge_been_search[index_closest_point] == true)
				continue;*/
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

			//add cutting planes to original mesh
			using Coord = double;
			using NN = uint32_t;
			using PPoint = std::array<Coord, 2>;
			std::vector<std::vector<PPoint>> polygon(1);
			/*Eigen::MatrixXd temp_Vis;
			temp_Vis.resize(3, 1);
			temp_Vis(0, 0) = Selected_Planes[i].origin.x();
			temp_Vis(1, 0) = Selected_Planes[i].origin.y();
			temp_Vis(2, 0) = Selected_Planes[i].origin.z();
			temp_Vis = rotMatrix.inverse() * temp_Vis;
			vis_cutting_planes[i].push_back(Vec3(temp_Vis(0, 0) + 50, temp_Vis(1, 0) + 50, temp_Vis(2, 0)));
			vis_cutting_planes[i].push_back(Vec3(temp_Vis(0, 0) + 50, temp_Vis(1, 0) - 50, temp_Vis(2, 0)));
			vis_cutting_planes[i].push_back(Vec3(temp_Vis(0, 0) - 50, temp_Vis(1, 0) - 50, temp_Vis(2, 0)));
			vis_cutting_planes[i].push_back(Vec3(temp_Vis(0, 0) - 50, temp_Vis(1, 0) + 50, temp_Vis(2, 0)));*/

			//delete triangles above the cutting plane
			removed_triangles_from_slicer.clear();
			removed_triangles_from_slicer.positions = slicer.positions;
			for (int j = 0; j < slicer.triangles.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (temp_slicer.positions[slicer.triangles[j][k]][2] - temp_slicer.positions[cutting_plane_polygons[0][0]][2] > 0.0000001) {
						removed_triangles_from_slicer.triangles.push_back(slicer.triangles[j]);
						slicer.triangles.erase(slicer.triangles.begin() + j);
						j--;
						break;
					}
				}
			}
			removed_triangles_from_slicer.save("temp_vis\\"+ mesh_target +"\\removed_triangles-" + to_string(i) + ".obj");

			Slicer_2 slicer_without_cutting_planes = slicer;
			for (int t = 0; t < cutting_plane_polygons.size(); t++)
			{
				//calculate the barycenter
				Eigen::MatrixXd temp_barycenter;
				temp_barycenter.resize(3, 1);
				double barycenter_x = 0, barycenter_y = 0, Area = 0;
				for (int j = 0; j < cutting_plane_polygons[t].size()-1; j++) {
					double x_1 = temp_slicer.positions[cutting_plane_polygons[t][j]][0];
					double y_1 = temp_slicer.positions[cutting_plane_polygons[t][j]][1];
					double x_2 = temp_slicer.positions[cutting_plane_polygons[t][j+1]][0];
					double y_2 = temp_slicer.positions[cutting_plane_polygons[t][j + 1]][1];
					barycenter_x += (x_1 + x_2)*(x_1*y_2 - x_2 * y_1);
					barycenter_y += (y_1 + y_2)*(x_1*y_2 - x_2 * y_1);
					Area += x_1*y_2 - x_2 * y_1;
				}
				Area /= 2;
				barycenter_x /= 6 * Area;
				barycenter_y /= 6 * Area;
				temp_barycenter(0, 0) = barycenter_x;
				temp_barycenter(1, 0) = barycenter_y;
				temp_barycenter(2, 0) = temp_slicer.positions[cutting_plane_polygons[t][0]][2];
				Eigen::Matrix3d rotMatrix_2 = Eigen::Quaterniond::FromTwoVectors(vector_after,vectorBefore).toRotationMatrix();
				temp_barycenter = rotMatrix_2.inverse() * temp_barycenter;
				barycenter_of_cutting_planes.push_back(Point_3(temp_barycenter(0,0), temp_barycenter(1, 0), temp_barycenter(2, 0)));
				
				map<int, int> map_index_faces;
				map_index_faces.clear();
				polygon[0].clear();
				//Current_Planes.push_back(Selected_Planes[i]);
				vector<Point_3> temp_vec;
				//Polygons_of_Current_Planes.push_back(temp_vec);
				for (int j = 0; j < cutting_plane_polygons[t].size(); j++) {
					map_index_faces.insert({ j,cutting_plane_polygons[t][j] });
					PPoint temp_point;
					temp_point[0] = temp_slicer.positions[cutting_plane_polygons[t][j]][0];
					temp_point[1] = temp_slicer.positions[cutting_plane_polygons[t][j]][1];
					polygon[0].push_back(temp_point);
					/*Polygons_of_Current_Planes[Polygons_of_Current_Planes.size() - 1].push_back(Point_3(slicer.positions[cutting_plane_polygons[t][j]][0],
						slicer.positions[cutting_plane_polygons[t][j]][1], slicer.positions[cutting_plane_polygons[t][j]][2]));*/
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
			slicer.save("temp_vis\\" + mesh_target + "\\(C)new_D_I-" + to_string(i) + ".obj");
			slicer_without_cutting_planes.save("temp_vis\\" + mesh_target + "\\new_D_I-" + to_string(i) + ".obj");
			if(percent_offset == 0 && i == Selected_Planes.size()-1)
				slicer.save("temp_vis\\" + mesh_target + "\\before_offset.obj");
		}
		//检查可达性
		index_last_useful_planes = 0;
		for (int i = 1; i < Selected_Planes.size(); i++)
			if (flag_useful_cutting_planes[i] == true)
				index_last_useful_planes = i;
		Mesh new_Di;
		Eigen::MatrixXd temp_V2;
		Eigen::MatrixXi temp_F2;
		igl::readOBJ("temp_vis\\" + mesh_target + "\\(C)new_D_I-" + to_string(index_last_useful_planes) + ".obj", temp_V2, temp_F2);
		igl::writeOFF("temp_vis\\" + mesh_target + "\\(C)new_D_I-" + to_string(index_last_useful_planes) + ".off", temp_V2, temp_F2);
		CGAL::IO::read_OFF("temp_vis\\" + mesh_target + "\\(C)new_D_I-" + to_string(index_last_useful_planes) + ".off", new_Di);
		std::vector<int> index_removed_Di_points;
		std::vector<Eigen::MatrixXd> vis_red_points_DS;
		std::vector<Eigen::MatrixXd> vis_red_points_DA;
		vis_red_points_DS.clear(); vis_red_points_DA.clear();
		int itr = percent_offset/step;    
		CGAL::Side_of_triangle_mesh<Mesh, K> inside(new_Di);
		for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
			Point_3 p(sampling_points_in_D_I[i].x(), sampling_points_in_D_I[i].y(), sampling_points_in_D_I[i].z());
			CGAL::Bounded_side res = inside(p);
			if (res != CGAL::ON_BOUNDED_SIDE)
				index_removed_Di_points.push_back(i);
		}
		Update_scalar_field_DS_DA(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
		Update_scalar_field_DS_DA_with_removed_points(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);  
		if (vis_red_points_DS.size() <= 1 && vis_red_points_DA.size() <= 1 && percent_offset >= step) //为了提高后续的搜索空间，强制要求至少offset一次    //底部的可达性检测还有问题，正常应该是==0时break; 13    //hand 100
			break;
		jud_is_force_offset = false;
		if (vis_red_points_DS.size() <= 1 && vis_red_points_DA.size() <= 1 && percent_offset < step) {
			jud_is_force_offset = true;
		}
		std::cout << percent_offset / step + 1 << " times of offset" << endl;
	}
	std::cout << "Initializing cutting planes over !!!" << endl;

	//visualize initial cutting planes
	/*ofstream ofile_2("temp_vis/initial cutting planes.obj");
	for (int i = 0; i < Selected_Planes.size(); i++) {
		Eigen::Vector3d vectorBefore(Selected_Planes[i].normal.x(), Selected_Planes[i].normal.y(), Selected_Planes[i].normal.z());
		Eigen::Vector3d vector_after(0,0,1);
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
		Eigen::MatrixXd temp_Vis;
		for (int j = 0; j < vis_cutting_planes[i].size(); j++) {
			temp_Vis.resize(3, 1);
			temp_Vis(0, 0) = vis_cutting_planes[i][j].m_x;
			temp_Vis(1, 0) = vis_cutting_planes[i][j].m_y;
			temp_Vis(2, 0) = vis_cutting_planes[i][j].m_z;
			temp_Vis = rotMatrix.inverse() * temp_Vis;
			vis_cutting_planes[i][j].m_x = temp_Vis(0, 0);
			vis_cutting_planes[i][j].m_y = temp_Vis(1, 0);
			vis_cutting_planes[i][j].m_z = temp_Vis(2, 0);
			ofile_2 << "v " << vis_cutting_planes[i][j].m_x << " " << vis_cutting_planes[i][j].m_y << " " << vis_cutting_planes[i][j].m_z << endl;
		}
		ofile_2 << "f " << to_string(i * 4+1+1) <<" " << to_string(i * 4+1) << " " << to_string(i * 4 + 2 + 1) << endl;
		ofile_2 << "f " << to_string(i * 4+3 + 1) << " " << to_string(i * 4 + 2+ 1) << " " << to_string(i * 4 + 1) << endl;
	}*/
	
	

}

void ReFab::Calculate_mean_curvature(string additive_interface,bool vis, double& mean_mean_principal_curvature, double& abs_mean_mean_principal_curvature)
{
	mean_mean_principal_curvature = 0;
	abs_mean_mean_principal_curvature = 0;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	using namespace Eigen;
	// Load a mesh in OFF format
	igl::read_triangle_mesh(additive_interface, V, F);

	// Alternative discrete mean curvature
	MatrixXd HN;
	SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	// Laplace-Beltrami of position
	HN = -Minv * (L * V);
	// Extract magnitude as mean curvature
	VectorXd H = HN.rowwise().norm();

	// Compute curvature directions via quadric fitting
	MatrixXd PD1, PD2;
	VectorXd PV1, PV2;
	igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
	// mean curvature
	H = 0.5 * (PV1 + PV2);
	for (int i = 0; i < H.rows(); i++) {
		mean_mean_principal_curvature += H(i, 0);
		abs_mean_mean_principal_curvature += abs(H(i, 0));
	}
	mean_mean_principal_curvature /= H.rows();
	//igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(V, F);

	//viewer.data().set_data(H);

	//// Average edge length for sizing
	//const double avg = igl::avg_edge_length(V, F);

	//// Draw a red segment parallel to the maximal curvature direction
	//const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
	//viewer.data().add_edges(V + PD1 * avg, V - PD1 * avg, red);

	//// Draw a blue segment parallel to the minimal curvature direction
	//viewer.data().add_edges(V + PD2 * avg, V - PD2 * avg, blue);

	//// Hide wireframe
	//viewer.data().show_lines = false;

	//viewer.launch();
}

void ReFab::Modify_Cutting_Planes()
{
	//提取每个切平面移除mesh的每个连通分量
	Eigen::MatrixXd temp_V;
	Eigen::MatrixXi temp_F;
	for (int i = 0; i < Selected_Planes.size(); i++) {
		igl::readOBJ("temp_vis\\removed_triangles-" + to_string(i) + ".obj", temp_V, temp_F);
		igl::writeOFF("temp_vis\\removed_triangles-" + to_string(i) + ".off", temp_V, temp_F);
		Mesh mesh;
		CGAL::IO::read_OFF("temp_vis\\removed_triangles-" + to_string(i) + ".off", mesh);
		std::vector<std::vector<int>> connected_components_triangles;
		connected_components_triangles.clear();
		vector<bool> flag_face_been_search(mesh.number_of_faces());
		for (int j = 0; j < mesh.number_of_faces(); j++)
			flag_face_been_search[j] = false;
		typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
		const double bound = std::cos(0.75 * CGAL_PI);
		std::vector<face_descriptor> cc;
		for (face_descriptor fd : faces(mesh))
		{
			cc.clear();
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
		Slicer_2 ori_slicer;
		ori_slicer.clear();
		ori_slicer.load("temp_vis\\new_D_I-" + to_string(i) + ".obj");
		for (int com = 0; com < connected_components_triangles.size(); com++) {
			Slicer_2 current_connected_component_mesh;
			current_connected_component_mesh.clear();
			current_connected_component_mesh.positions.resize(temp_V.rows());
			for (int j = 0; j < temp_V.rows(); j++) {
				current_connected_component_mesh.positions[j][0] = temp_V(j, 0);
				current_connected_component_mesh.positions[j][1] = temp_V(j, 1);
				current_connected_component_mesh.positions[j][2] = temp_V(j, 2);
			}
			current_connected_component_mesh.triangles.resize(connected_components_triangles[com].size());
			for (int j = 0; j < connected_components_triangles[com].size(); j++) {
				current_connected_component_mesh.triangles[j][0] = temp_F(connected_components_triangles[com][j],0);
				current_connected_component_mesh.triangles[j][1] = temp_F(connected_components_triangles[com][j],1);
				current_connected_component_mesh.triangles[j][2] = temp_F(connected_components_triangles[com][j],2);
			}
			double offset;
			offset = 2;
			if(current_connected_component_mesh.triangles.size()<300)   //暂时将小于300个三角形的连通分量的切平面偏移量设为0.2，需要修改
				offset = 0.02;
			Plane current_new_plane;    //需要修改切平面的生成方式，且需要添加必须是contour的约束
			current_new_plane.normal = Point_3(Selected_Planes[i].normal[0], Selected_Planes[i].normal[1], Selected_Planes[i].normal[2]);
			current_new_plane.origin = Point_3(Selected_Planes[i].origin.x() + offset * Selected_Planes[i].normal[0], Selected_Planes[i].origin.y() + offset * Selected_Planes[i].normal[1], 
				Selected_Planes[i].origin.z() + offset * Selected_Planes[i].normal[2]);
			current_connected_component_mesh.normal[0] = current_new_plane.normal.x();
			current_connected_component_mesh.normal[1] = current_new_plane.normal.y();
			current_connected_component_mesh.normal[2] = current_new_plane.normal.z();
			current_connected_component_mesh.origin[0] = current_new_plane.origin.x();
			current_connected_component_mesh.origin[1] = current_new_plane.origin.y();
			current_connected_component_mesh.origin[2] = current_new_plane.origin.z();
			current_connected_component_mesh.cut();
			//current_connected_component_mesh.save("temp_vis\\11-.obj");
			ori_slicer.positions = current_connected_component_mesh.positions;
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(current_connected_component_mesh.normal[0], current_connected_component_mesh.normal[1], current_connected_component_mesh.normal[2]);
			Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			Slicer_2 temp_slicer = current_connected_component_mesh;
			for (int j = 0; j < temp_slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = temp_slicer.positions[j][0];
				temp_V(1, 0) = temp_slicer.positions[j][1];
				temp_V(2, 0) = temp_slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				temp_slicer.positions[j][0] = temp_V(0, 0);
				temp_slicer.positions[j][1] = temp_V(1, 0);
				temp_slicer.positions[j][2] = temp_V(2, 0);
			}
			Eigen::MatrixXd temp_origin;
			temp_origin.resize(3, 1);
			temp_origin(0, 0) = current_connected_component_mesh.origin[0];
			temp_origin(1, 0) = current_connected_component_mesh.origin[1];
			temp_origin(2, 0) = current_connected_component_mesh.origin[2];
			temp_origin = rotMatrix.inverse() * temp_origin;
			Slicer_2 removed_triangles_from_slicer;
			removed_triangles_from_slicer.clear();
			removed_triangles_from_slicer.positions = current_connected_component_mesh.positions;
			for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (temp_slicer.positions[current_connected_component_mesh.triangles[j][k]][2] - temp_origin(2, 0) > 0.0000001) {
						removed_triangles_from_slicer.triangles.push_back(current_connected_component_mesh.triangles[j]);
						current_connected_component_mesh.triangles.erase(current_connected_component_mesh.triangles.begin() + j);
						j--;
						break;
					}
				}
			}
			temp_V.resize(current_connected_component_mesh.positions.size(), 3);
			for (int j = 0; j < current_connected_component_mesh.positions.size(); j++) {
				temp_V(j, 0) = current_connected_component_mesh.positions[j][0];
				temp_V(j, 1) = current_connected_component_mesh.positions[j][1];
				temp_V(j, 2) = current_connected_component_mesh.positions[j][2];
			}

			//将当前连通分量保留的mesh加入到原始切平面中
			Add_triangles_of_cutting_planes(current_connected_component_mesh);
			for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
				ori_slicer.triangles.push_back(current_connected_component_mesh.triangles[j]);
			}
			removed_triangles_from_slicer.save("temp_vis\\removed_triangles-" + to_string(com) + "itr-" +to_string(1) + ".obj");
			ori_slicer.save("temp_vis\\temp_new_Di.obj");

			//计算此时的可达性情况
			Mesh new_Di;
			Eigen::MatrixXd temp_V2;
			Eigen::MatrixXi temp_F2;
			igl::readOBJ("temp_vis\\temp_new_Di.obj", temp_V2, temp_F2);
			igl::writeOFF("temp_vis\\temp_new_Di.off", temp_V2, temp_F2);
			CGAL::IO::read_OFF("temp_vis\\temp_new_Di.off", new_Di);
			std::vector<int> index_removed_Di_points;
			std::vector<Eigen::MatrixXd> vis_red_points_DS;
			std::vector<Eigen::MatrixXd> vis_red_points_DA;
			vis_red_points_DS.clear(); vis_red_points_DA.clear();
			int itr = Selected_Planes.size() - i;
			CGAL::Side_of_triangle_mesh<Mesh, K> inside(new_Di);
			for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
				Point_3 p(sampling_points_in_D_I[i].x(), sampling_points_in_D_I[i].y(), sampling_points_in_D_I[i].z());
				CGAL::Bounded_side res = inside(p);
				if (res != CGAL::ON_BOUNDED_SIDE)
					index_removed_Di_points.push_back(i);
			}
			Update_scalar_field_DS_DA(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
			Update_scalar_field_DS_DA_with_removed_points(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
			//计算此时的平凸约束（先不考虑顶点位置导致的凸变凹）
			std::vector<Point_3> additive_interface_points_inside_Di;
			for (int i = 0; i < additive_interface_points.size(); i++) {
				for (int j = 0; j < additive_interface_points[i].size(); j++) {
					Point_3 p = additive_interface_points[i][j];
					CGAL::Bounded_side res = inside(p);
					if (res == CGAL::ON_BOUNDED_SIDE)
						additive_interface_points_inside_Di.push_back(p);
				}
			}
			if (additive_interface_points_inside_Di.size() == 0)
				cout << "Current cutting planes are satisfy the flatten cosntraint √√√" << endl;
			else
				cout << "Current cutting planes aren't satisfy the flatten cosntraint ×××" << endl;
		}
		ori_slicer.save("temp_vis\\(C)new_D_I-" + to_string(i) + ".obj");

		//更新下一个切平面
		if (i < Selected_Planes.size() - 1) {
			while (flag_use_selected_planes[i + 1] == false) {
				i++;
				if (i == Selected_Planes.size() - 1)
					break;
			}
			ori_slicer.normal[0] = Selected_Planes[i+1].normal[0];
			ori_slicer.normal[1] = Selected_Planes[i+1].normal[1];
			ori_slicer.normal[2] = Selected_Planes[i+1].normal[2];
			ori_slicer.origin[0] = Selected_Planes[i+1].origin.x();
			ori_slicer.origin[1] = Selected_Planes[i+1].origin.y();
			ori_slicer.origin[2] = Selected_Planes[i+1].origin.z();
			ori_slicer.cut();
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(ori_slicer.normal[0], ori_slicer.normal[1], ori_slicer.normal[2]);
			Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int j = 0; j < ori_slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = ori_slicer.positions[j][0];
				temp_V(1, 0) = ori_slicer.positions[j][1];
				temp_V(2, 0) = ori_slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				ori_slicer.positions[j][0] = temp_V(0, 0);
				ori_slicer.positions[j][1] = temp_V(1, 0);
				ori_slicer.positions[j][2] = temp_V(2, 0);
			}
			Eigen::MatrixXd temp_origin;
			temp_origin.resize(3, 1);
			temp_origin(0, 0) = ori_slicer.origin[0];
			temp_origin(1, 0) = ori_slicer.origin[1];
			temp_origin(2, 0) = ori_slicer.origin[2];
			temp_origin = rotMatrix.inverse() * temp_origin;
			Slicer_2 removed_triangles_from_slicer, remain_triangles_form_slicer;
			removed_triangles_from_slicer.clear(); remain_triangles_form_slicer.clear();
			for (int j = 0; j < ori_slicer.triangles.size(); j++) {
				bool jud_removed = false;
				for (int k = 0; k < 3; k++) {
					if (ori_slicer.positions[ori_slicer.triangles[j][k]][2] - temp_origin(2, 0) > 0.0000001) {
						removed_triangles_from_slicer.triangles.push_back(ori_slicer.triangles[j]);
						jud_removed = true;
						break;
					}
				}
				if (jud_removed == false)
					remain_triangles_form_slicer.triangles.push_back(ori_slicer.triangles[j]);
			}
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after,vectorBefore).toRotationMatrix();
			for (int j = 0; j < ori_slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = ori_slicer.positions[j][0];
				temp_V(1, 0) = ori_slicer.positions[j][1];
				temp_V(2, 0) = ori_slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				ori_slicer.positions[j][0] = temp_V(0, 0);
				ori_slicer.positions[j][1] = temp_V(1, 0);
				ori_slicer.positions[j][2] = temp_V(2, 0);
			}
			removed_triangles_from_slicer.positions = ori_slicer.positions;
			remain_triangles_form_slicer.positions = ori_slicer.positions;
			removed_triangles_from_slicer.save("temp_vis\\removed_triangles-" + to_string(i + 1) + ".obj");
			remain_triangles_form_slicer.save("temp_vis\\new_D_I-" + to_string(i + 1) + ".obj");
		}
	}
}

bool ReFab::Check_is_one_component(Mesh new_Di)
{
	vector<bool> flag_face_been_search(new_Di.number_of_faces());
	for (int j = 0; j < new_Di.number_of_faces(); j++)
		flag_face_been_search[j] = false;
	std::vector<face_descriptor> cc;
	int cont_components = 0;
	for (face_descriptor fd : faces(new_Di))
	{
		cc.clear();
		if (flag_face_been_search[fd.idx()] == true)
			continue;
		PMP::connected_component(fd,
			new_Di,
			std::back_inserter(cc));
		std::vector<int> temp_points;
		for (int j = 0; j < cc.size(); j++) {
			flag_face_been_search[cc[j].idx()] = true;
		}
		cont_components++;
		if (cont_components >= 2)
			return false;
	}
	return true;
}


double ReFab::Calculate_distance_plane_and_point(Point_3 normal, Point_3 origin, Point_3 p)
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

bool ReFab::Modify_Cutting_Planes_2(vector<Plane> parent_planes,double& new_volume, vector<Plane>& best_modifed_cutting_planes, double& max_modified_volume, double& last_max_volume,vector<Point_3>& new_current_removed_points,bool flag_final, bool& is_max_volume)
{
	//提取每个切平面移除mesh的每个连通分量
	ofstream save_final_cutting_planes("MCTS_temp\\" + mesh_target + "\\final_cutting_planes.txt");
	Eigen::MatrixXd temp_V;
	Eigen::MatrixXi temp_F;
	int cont_planes = 0;
	int cont_com = 0;
	new_volume = 0;
	vector<Plane> save_current_modified_cutting_planes;
	save_current_modified_cutting_planes.clear();
	clock_t start_time = clock();
	double sum_time = 0;
	double sum_time2 = 0;
	double sum_time3 = 0;
	double sum_time4 = 0;
	double sum_file = 0;
	clock_t sf, ef;
	Slicer_2 ori_slicer, removed_slicer;
	ori_slicer.clear();
	int cont_final_cutting_planes = 0;
	vector<Plane> real_final_cutting_planes;
	Slicer_2 save_tiny_triangles;
	real_final_cutting_planes.clear();
	for (int i = 0; i < Selected_Planes.size(); i++) {
		removed_slicer = saved_removed_triangles;
		bool* flag_remain_removed_triangles = new bool [saved_removed_triangles.triangles.size()];
		for(int j= 0;j< saved_removed_triangles.triangles.size();j++)
			flag_remain_removed_triangles[j] = true;
		ori_slicer = saved_new_D_I;
		temp_V.resize(saved_removed_triangles.positions.size(),3);
		temp_F.resize(saved_removed_triangles.triangles.size(),3);
		for (int j = 0; j < saved_removed_triangles.positions.size(); j++) {
			temp_V(j, 0) = saved_removed_triangles.positions[j][0];
			temp_V(j, 1) = saved_removed_triangles.positions[j][1];
			temp_V(j, 2) = saved_removed_triangles.positions[j][2];
		}
		for (int j = 0; j < saved_removed_triangles.triangles.size(); j++) {
			temp_F(j,0) = saved_removed_triangles.triangles[j][0];
			temp_F(j,1) = saved_removed_triangles.triangles[j][1];
			temp_F(j,2) = saved_removed_triangles.triangles[j][2];
		}
			
		/*igl::readOBJ("MCTS_temp\\" + mesh_target + "\\removed_triangles-" + to_string(i) + ".obj", temp_V, temp_F);
		Mesh mesh;
		CGAL::IO::read_OBJ("MCTS_temp\\" + mesh_target + "\\removed_triangles-" + to_string(i) + ".obj", mesh);*/
		/*Slicer_2 ori_slicer;
		ori_slicer.clear();
		ori_slicer.load("MCTS_temp\\" + mesh_target + "\\new_D_I-" + to_string(i) + ".obj");*/

		
		/*vector<bool> flag_face_been_search(mesh.number_of_faces());
		for (int j = 0; j < mesh.number_of_faces(); j++)
			flag_face_been_search[j] = false;
		typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
		const double bound = std::cos(0.75 * CGAL_PI);
		std::vector<face_descriptor> cc;*/
		
		//需要加速
				/*for (face_descriptor fd : faces(mesh))
				{
					cc.clear();
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
				}*/
		clock_t st = clock();
		std::vector<std::vector<int>> connected_components_triangles;
		connected_components_triangles.clear();
		Eigen::VectorXi components;
		igl::facet_components(temp_F, components);
		int num_components = 0;
		for (int i = 0; i < components.rows(); ++i)
			num_components = std::max(num_components, components(i));
		connected_components_triangles.resize(num_components+1);
		for (int i = 0; i < components.rows(); ++i) {
			connected_components_triangles[components(i)].push_back(i);
		}
		clock_t et = clock();
		sum_time += (double)(et - st) / CLOCKS_PER_SEC;  

		clock_t st2 = clock();
		
		sf = clock();
		
		ef = clock();
		sum_file += (double)(ef - sf) / CLOCKS_PER_SEC;
		clock_t et2 = clock();
		sum_time2 += (double)(et2 - st2) / CLOCKS_PER_SEC;  
		for (int com = 0; com < connected_components_triangles.size(); com++) {
			Slicer_2 current_connected_component_mesh;
			current_connected_component_mesh.clear();
			current_connected_component_mesh.positions.resize(temp_V.rows());
			for (int j = 0; j < temp_V.rows(); j++) {
				current_connected_component_mesh.positions[j][0] = temp_V(j, 0);
				current_connected_component_mesh.positions[j][1] = temp_V(j, 1);
				current_connected_component_mesh.positions[j][2] = temp_V(j, 2);
			}
			current_connected_component_mesh.triangles.resize(connected_components_triangles[com].size());
			for (int j = 0; j < connected_components_triangles[com].size(); j++) {
				current_connected_component_mesh.triangles[j][0] = temp_F(connected_components_triangles[com][j], 0);
				current_connected_component_mesh.triangles[j][1] = temp_F(connected_components_triangles[com][j], 1);
				current_connected_component_mesh.triangles[j][2] = temp_F(connected_components_triangles[com][j], 2);
				flag_remain_removed_triangles[connected_components_triangles[com][j]] = false;
			}
			Plane current_new_plane;  
			bool flag_last_plane_determine;
			if (cont_planes < parent_planes.size()) {
				current_new_plane.normal = parent_planes[cont_planes].normal;
				current_new_plane.origin = parent_planes[cont_planes].origin;
				flag_last_plane_determine = true;
			}
			else {
				//std::vector<Point_3> sample_normals = OrientationSamplePoints_normal();
				//int random_index = rand() % sample_normals.size();
				//current_new_plane.normal = sample_normals[random_index];
				current_new_plane.normal = Selected_Planes[i].normal;
				Point_3 origin_point(barycenter_of_cutting_planes[cont_planes].x(), barycenter_of_cutting_planes[cont_planes].y(), barycenter_of_cutting_planes[cont_planes].z());
				double max_offset_distance;
				if (index_of_additive_interface[cont_planes].size() == 0) {  //求最大距离
					double max_distance = -100000000;
					if (i == 0) {
						for (int t = 0; t < all_vertices_in_components[cont_planes].size(); t++) {
							double distance = Calculate_distance_plane_and_point(current_new_plane.normal, origin_point, all_vertices_in_components[cont_planes][t]);
							if (distance > max_distance) {
								max_distance = distance;
							}
						}
					}
					else {   //更新点
						vector<Point_3> new_current_points;
						new_current_points.clear();
						for (int t = 0; t < temp_F.rows(); t++)
							for (int j = 0; j < 3; j++) 
								new_current_points.push_back(Point_3(temp_V(temp_F(t, j), 0), temp_V(temp_F(t, j), 1), temp_V(temp_F(t, j), 2)));
						for (int t = 0; t < new_current_points.size(); t++) {
							double distance = Calculate_distance_plane_and_point(current_new_plane.normal, origin_point, new_current_points[t]);
							if (distance > max_distance) {
								max_distance = distance;
							}
						}
						if(flag_last_plane_determine == true)
							new_current_removed_points = new_current_points;
					}
					max_offset_distance = max_distance;
				}
				else {  //求最小距离
					double min_distance = 100000000;
					for (int t = 0; t < index_of_additive_interface[cont_planes].size(); t++) {
						for (int j = 0; j < additive_interface_points[index_of_additive_interface[cont_planes][t]].size(); j++) {
							double distance = Calculate_distance_plane_and_point(current_new_plane.normal, origin_point, additive_interface_points[index_of_additive_interface[cont_planes][t]][j]);
							if (distance < min_distance) {
								min_distance = distance;
							}
						}
					}
					max_offset_distance = min_distance;
				}
				float random_offset = 0.0f + static_cast<float>(rand()) / RAND_MAX * max_offset_distance;   //伪随机数
				//float random_offset = max_offset_distance - 0.1;
				
				current_new_plane.origin = Point_3(barycenter_of_cutting_planes[cont_planes].x() + current_new_plane.normal.x() * random_offset,
					barycenter_of_cutting_planes[cont_planes].y() + current_new_plane.normal.y() * random_offset, barycenter_of_cutting_planes[cont_planes].z() + current_new_plane.normal.z() * random_offset);
				flag_last_plane_determine = false;
			}
			if (best_modifed_cutting_planes.size() != 0) {
				current_new_plane.normal = best_modifed_cutting_planes[cont_com].normal;
				current_new_plane.origin = best_modifed_cutting_planes[cont_com].origin;
			}
			else
			{
				save_current_modified_cutting_planes.push_back(current_new_plane);
			}
			cont_com++;
			
			int num_triangles_before_cutting = current_connected_component_mesh.positions.size();
			bool flag_tiny_component = false;
			current_connected_component_mesh.normal[0] = current_new_plane.normal.x();
			current_connected_component_mesh.normal[1] = current_new_plane.normal.y();
			current_connected_component_mesh.normal[2] = current_new_plane.normal.z();
			current_connected_component_mesh.origin[0] = current_new_plane.origin.x();
			current_connected_component_mesh.origin[1] = current_new_plane.origin.y();
			current_connected_component_mesh.origin[2] = current_new_plane.origin.z();
			current_connected_component_mesh.cut();
			ori_slicer.positions = current_connected_component_mesh.positions;
			if (num_triangles_before_cutting == current_connected_component_mesh.positions.size())
				flag_tiny_component = true;
			if (flag_tiny_component == true) {
				for(int j = 0; j < current_connected_component_mesh.triangles.size(); j++)
					save_tiny_triangles.triangles.push_back(current_connected_component_mesh.triangles[j]);
			}		
			if (flag_tiny_component == false) {
				Eigen::Vector3d vectorBefore(0, 0, 1);
				Eigen::Vector3d vector_after(current_connected_component_mesh.normal[0], current_connected_component_mesh.normal[1], current_connected_component_mesh.normal[2]);
				Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
				Slicer_2 temp_slicer = current_connected_component_mesh;
				for (int j = 0; j < temp_slicer.positions.size(); j++) {
					Eigen::MatrixXd temp_V;
					temp_V.resize(3, 1);
					temp_V(0, 0) = temp_slicer.positions[j][0];
					temp_V(1, 0) = temp_slicer.positions[j][1];
					temp_V(2, 0) = temp_slicer.positions[j][2];
					temp_V = rotMatrix.inverse() * temp_V;
					temp_slicer.positions[j][0] = temp_V(0, 0);
					temp_slicer.positions[j][1] = temp_V(1, 0);
					temp_slicer.positions[j][2] = temp_V(2, 0);
				}
				
				Eigen::MatrixXd temp_origin;
				temp_origin.resize(3, 1);
				temp_origin(0, 0) = current_connected_component_mesh.origin[0];
				temp_origin(1, 0) = current_connected_component_mesh.origin[1];
				temp_origin(2, 0) = current_connected_component_mesh.origin[2];
				temp_origin = rotMatrix.inverse() * temp_origin;
				Slicer_2 removed_triangles_from_slicer;
				removed_triangles_from_slicer.clear();
				removed_triangles_from_slicer.positions = current_connected_component_mesh.positions;


				//clock_t st4 = clock();
				//0.08s需要加速
						for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
							for (int k = 0; k < 3; k++) {
								if (temp_slicer.positions[current_connected_component_mesh.triangles[j][k]][2] - temp_origin(2, 0) > 0.0000001) {
									removed_triangles_from_slicer.triangles.push_back(current_connected_component_mesh.triangles[j]);
									current_connected_component_mesh.triangles.erase(current_connected_component_mesh.triangles.begin() + j);
									j--;
									break;
								}
							}
						}
				//clock_t et4 = clock();
				//sum_time4 += (double)(et4 - st4) / CLOCKS_PER_SEC;  
				double min_point_of_removed_triangles = 100000000;
				int index_of_min_triangles = -1;
				int index_of_min_triangles_2 = -1;
				for (int j = 0; j < removed_triangles_from_slicer.triangles.size(); j++) {
					for (int k = 0; k < 3; k++) {
						if (removed_triangles_from_slicer.positions[removed_triangles_from_slicer.triangles[j][k]][2] < min_point_of_removed_triangles) {
							min_point_of_removed_triangles = removed_triangles_from_slicer.positions[removed_triangles_from_slicer.triangles[j][k]][2];
							index_of_min_triangles = j;
							index_of_min_triangles_2 = k;
						}
					}
				}
				temp_V.resize(current_connected_component_mesh.positions.size(), 3);
				for (int j = 0; j < current_connected_component_mesh.positions.size(); j++) {
					temp_V(j, 0) = current_connected_component_mesh.positions[j][0];
					temp_V(j, 1) = current_connected_component_mesh.positions[j][1];
					temp_V(j, 2) = current_connected_component_mesh.positions[j][2];
				}
				//removed_triangles_from_slicer.save("MCTS_temp\\removed_triangles_ttttmp.obj");
				if (flag_final == true) { //存放当前cutting plane修改后切割的结果
					//removed_triangles_from_slicer.save("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(i) +"-"+to_string(com)+ ".obj");
					removed_triangles_from_slicer.save("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(cont_final_cutting_planes) + ".obj");
					final_cutting_planes.push_back(Plane(Point_3(current_new_plane.normal.x(), current_new_plane.normal.y(), current_new_plane.normal.z()),
						Point_3(current_new_plane.origin.x(), current_new_plane.origin.y(), current_new_plane.origin.z())));
					points_in_cutting_planes.push_back(Point_3(removed_triangles_from_slicer.positions[removed_triangles_from_slicer.triangles[index_of_min_triangles][index_of_min_triangles_2]][0], 
						removed_triangles_from_slicer.positions[removed_triangles_from_slicer.triangles[index_of_min_triangles][index_of_min_triangles_2]][1], 
						removed_triangles_from_slicer.positions[removed_triangles_from_slicer.triangles[index_of_min_triangles][index_of_min_triangles_2]][2]));
					save_final_cutting_planes << current_new_plane.normal.x() << " " << current_new_plane.normal.y() << " " << current_new_plane.normal.z() << " " << current_new_plane.origin.x() << " " << current_new_plane.origin.y() << " " << current_new_plane.origin.z() << endl;
					save_final_cutting_planes << points_in_cutting_planes[points_in_cutting_planes.size() - 1].x() << " " << points_in_cutting_planes[points_in_cutting_planes.size() - 1].y() << " " << points_in_cutting_planes[points_in_cutting_planes.size() - 1].z() << endl;
					real_final_cutting_planes.push_back(best_modifed_cutting_planes[cont_com - 1]);
				
					Slicer_2 ori_slicer_2 = ori_slicer;
					for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
						ori_slicer_2.triangles.push_back(current_connected_component_mesh.triangles[j]);
					}
					for (int j = 0; j < saved_removed_triangles.triangles.size(); j++) {
						if (flag_remain_removed_triangles[j] == true) {
							ori_slicer_2.triangles.push_back(saved_removed_triangles.triangles[j]);
						}
					}
					/*for(int j =0;j< save_tiny_triangles.triangles.size();j++)
						ori_slicer_2.triangles.push_back(save_tiny_triangles.triangles[j]);*/
					save_tiny_triangles.triangles.clear();
					/*ori_slicer_2.normal[0] = current_new_plane.normal.x();
					ori_slicer_2.normal[1] = current_new_plane.normal.y();
					ori_slicer_2.normal[2] = current_new_plane.normal.z();
					ori_slicer_2.origin[0] = current_new_plane.origin.x();
					ori_slicer_2.origin[1] = current_new_plane.origin.y();
					ori_slicer_2.origin[2] = current_new_plane.origin.z();
					Add_triangles_of_cutting_planes(ori_slicer_2);*/
					ori_slicer_2.save("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(cont_final_cutting_planes) + ".obj");
				}
				cont_final_cutting_planes++;

				//将当前连通分量保留的mesh加入到原始切平面中
				if (!Add_triangles_of_cutting_planes(current_connected_component_mesh))
					return false;
				
				//if (flag_final == true) { //存放当前cutting plane修改后切割的结果
				//	ori_slicer.save("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(cont_planes) + ".obj");
				//}
			}
			
			for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
				ori_slicer.triangles.push_back(current_connected_component_mesh.triangles[j]);
			}
			
			if (i == Selected_Planes.size() - 1 && com == connected_components_triangles.size() - 1) {  //debug用，存放当前mesh
				ori_slicer.save("MCTS_temp\\" + mesh_target + "\\temp_new_Di.obj");
			}
				
			//if (flag_final == true && flag_tiny_component == false) { //存放当前cutting plane修改后切割的结果
			//	ori_slicer.save("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(cont_planes) + "aa.obj");
			//}

			//当修改后的cutting plane切割出的连通分量多于初始化cutting plane切割出的连通分量，则修改后得到的连通分量用的cutting plane数量向下兼容；
			//反之，则用的cutting plane数量向上兼容,这种情况还没有考虑怎么处理。。。。。。
			if(cont_planes < cont_components_of_each_planes[i] && cont_planes < cont_components_of_each_planes[cont_components_of_each_planes.size() - 1] - 1)
				cont_planes++;
		}
		
		
		//ori_slicer.save("MCTS_temp\\(C)new_D_I-" + to_string(i) + ".obj");
		//if (flag_final == true) { //存放当前cutting plane修改后切割的结果
		//	ori_slicer.save("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(i) + ".obj");
		//}
		
		//更新下一个切平面
		if (i < Selected_Planes.size() - 1) {
			/*while (flag_use_selected_planes[i + 1] == false) {
				i++;
				if (i == Selected_Planes.size() - 1)
					break;
			}*/
			clock_t st3 = clock();
			ori_slicer.normal[0] = Selected_Planes[i + 1].normal[0];
			ori_slicer.normal[1] = Selected_Planes[i + 1].normal[1];
			ori_slicer.normal[2] = Selected_Planes[i + 1].normal[2];
			ori_slicer.origin[0] = Selected_Planes[i + 1].origin.x();
			ori_slicer.origin[1] = Selected_Planes[i + 1].origin.y();
			ori_slicer.origin[2] = Selected_Planes[i + 1].origin.z();
			ori_slicer.cut();
			//ori_slicer.save("MCTS_temp\\aaaatmp.obj");
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vector_after(ori_slicer.normal[0], ori_slicer.normal[1], ori_slicer.normal[2]);
			Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
			for (int j = 0; j < ori_slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = ori_slicer.positions[j][0];
				temp_V(1, 0) = ori_slicer.positions[j][1];
				temp_V(2, 0) = ori_slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				ori_slicer.positions[j][0] = temp_V(0, 0);
				ori_slicer.positions[j][1] = temp_V(1, 0);
				ori_slicer.positions[j][2] = temp_V(2, 0);
			}
			Eigen::MatrixXd temp_origin;
			temp_origin.resize(3, 1);
			temp_origin(0, 0) = ori_slicer.origin[0];
			temp_origin(1, 0) = ori_slicer.origin[1];
			temp_origin(2, 0) = ori_slicer.origin[2];
			temp_origin = rotMatrix.inverse() * temp_origin;
			Slicer_2 removed_triangles_from_slicer, remain_triangles_form_slicer;
			removed_triangles_from_slicer.clear(); remain_triangles_form_slicer.clear();
			for (int j = 0; j < ori_slicer.triangles.size(); j++) {
				bool jud_removed = false;
				for (int k = 0; k < 3; k++) {
					if (ori_slicer.positions[ori_slicer.triangles[j][k]][2] - temp_origin(2, 0) > 0.0000001) {
						removed_triangles_from_slicer.triangles.push_back(ori_slicer.triangles[j]);
						jud_removed = true;
						break;
					}
				}
				if (jud_removed == false)
					remain_triangles_form_slicer.triangles.push_back(ori_slicer.triangles[j]);
			}
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vectorBefore).toRotationMatrix();
			for (int j = 0; j < ori_slicer.positions.size(); j++) {
				Eigen::MatrixXd temp_V;
				temp_V.resize(3, 1);
				temp_V(0, 0) = ori_slicer.positions[j][0];
				temp_V(1, 0) = ori_slicer.positions[j][1];
				temp_V(2, 0) = ori_slicer.positions[j][2];
				temp_V = rotMatrix.inverse() * temp_V;
				ori_slicer.positions[j][0] = temp_V(0, 0);
				ori_slicer.positions[j][1] = temp_V(1, 0);
				ori_slicer.positions[j][2] = temp_V(2, 0);
			}
			removed_triangles_from_slicer.positions = ori_slicer.positions;
			remain_triangles_form_slicer.positions = ori_slicer.positions;
			
			//0.25s需要加速
					sf = clock();
					saved_removed_triangles = removed_triangles_from_slicer;
					saved_new_D_I = remain_triangles_form_slicer;

					//removed_triangles_from_slicer.save("MCTS_temp\\" + mesh_target + "\\removed_triangles-" + to_string(i + 1) + ".obj");
					//remain_triangles_form_slicer.save("MCTS_temp\\" + mesh_target + "\\new_D_I-" + to_string(i + 1) + ".obj");
					ef = clock();
					sum_file += (double)(ef - sf) / CLOCKS_PER_SEC;
			clock_t et3 = clock();  
			sum_time3 += (double)(et3 - st3) / CLOCKS_PER_SEC; 
		}
	}
	//std::cout << "AAD time: " << sum_time4 << "s" << std::endl;
	//std::cout << "AAC time: " << sum_time3 << "s" << std::endl;
	//std::cout << "AAB time: " << sum_time2 << "s" << std::endl;
	//std::cout << "AAA time: " << sum_time << "s" << std::endl;
	clock_t end_time = clock();
	//std::cout << "AA time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;  
	//计算此时的可达性情况
	clock_t start_time2, end_time2;
	start_time2 = clock();
	
	clock_t start_time5 = clock();
	Mesh new_Di;
	Eigen::MatrixXd temp_V2;
	Eigen::MatrixXi temp_F2;
	temp_V2.resize(ori_slicer.positions.size(), 3);
	temp_F2.resize(ori_slicer.triangles.size(), 3);
	for (int j = 0; j < ori_slicer.positions.size(); j++) {
		temp_V2(j, 0) = ori_slicer.positions[j][0];
		temp_V2(j, 1) = ori_slicer.positions[j][1];
		temp_V2(j, 2) = ori_slicer.positions[j][2];
	}
	for (int j = 0; j < ori_slicer.triangles.size(); j++) {
		temp_F2(j, 0) = ori_slicer.triangles[j][0];
		temp_F2(j, 1) = ori_slicer.triangles[j][1];
		temp_F2(j, 2) = ori_slicer.triangles[j][2];
	}
	//igl::readOBJ("MCTS_temp\\" + mesh_target + "\\temp_new_Di.obj", temp_V2, temp_F2);
	//igl::writeOFF("MCTS_temp\\" + mesh_target + "\\temp_new_Di.off", temp_V2, temp_F2);
	//CGAL::IO::read_OBJ("MCTS_temp\\" + mesh_target + "\\temp_new_Di.obj", new_Di);
	/*if (!Check_is_one_component(new_Di))  //先注释，遇到问题再说
		return false;*/
	std::vector<int> index_removed_Di_points;
	std::vector<Eigen::MatrixXd> vis_red_points_DS;
	std::vector<Eigen::MatrixXd> vis_red_points_DA;
	vis_red_points_DS.clear(); vis_red_points_DA.clear();
	clock_t end_time5 = clock();
	//std::cout << "ABA time: " << (double)(end_time5 - start_time5) / CLOCKS_PER_SEC << "s" << std::endl;

	for (int i = 0; i < temp_V2.rows(); ++i) {
		Point_3 p(temp_V2(i, 0), temp_V2(i, 1), temp_V2(i, 2));
		new_Di.add_vertex(p);
	}
	for (int i = 0; i < temp_F2.rows(); ++i) {
		new_Di.add_face(
			Mesh::Vertex_index(temp_F2(i, 0)),
			Mesh::Vertex_index(temp_F2(i, 1)),
			Mesh::Vertex_index(temp_F2(i, 2))
		);
	}

	int itr = -1;
	CGAL::Side_of_triangle_mesh<Mesh, K> inside(new_Di);
	for (int i = 0; i < sampling_points_in_D_I.size(); i++) {
		Point_3 p(sampling_points_in_D_I[i].x(), sampling_points_in_D_I[i].y(), sampling_points_in_D_I[i].z());
		CGAL::Bounded_side res = inside(p);
		if (res != CGAL::ON_BOUNDED_SIDE)
			index_removed_Di_points.push_back(i);
	}
	clock_t end_time4 = clock();
	//std::cout << "ABB time: " << (double)(end_time4 - start_time4) / CLOCKS_PER_SEC << "s" << std::endl;

	clock_t start_time6 = clock();
	Update_scalar_field_DS_DA_2(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);
	clock_t end_time6 = clock();
	//std::cout << "ABC time: " << (double)(end_time6 - start_time6) / CLOCKS_PER_SEC << "s" << std::endl;
	clock_t start_time3 = clock();
	Update_scalar_field_DS_DA_with_removed_points_2(index_removed_Di_points, vis_red_points_DS, vis_red_points_DA, itr);  //itr取1可可视化不可达点
	clock_t end_time3 = clock();
	//std::cout << "ABD time: " << (double)(end_time3 - start_time3) / CLOCKS_PER_SEC << "s" << std::endl;
	end_time2 = clock();
	//std::cout << "AB time: " << (double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;

	if (vis_red_points_DS.size() > 1 || vis_red_points_DA.size() > 1)  //正常应该是＞0，目前底部有bug   暂时注释，测试
		return false;

	start_time2 = clock();
	//计算此时的平凸约束（先不考虑顶点位置导致的凸变凹）  //*********由于offset的设置，这个检查似乎不需要了
	std::vector<Point_3> additive_interface_points_inside_Di;
	for (int i = 0; i < additive_interface_points.size(); i++) {
		for (int j = 0; j < additive_interface_points[i].size(); j++) {
			Point_3 p = additive_interface_points[i][j];
			CGAL::Bounded_side res = inside(p);
			if (res == CGAL::ON_BOUNDED_SIDE || res == CGAL::ON_BOUNDARY)
				additive_interface_points_inside_Di.push_back(p);
		}
	}
	if (additive_interface_points_inside_Di.size() != 0)
		return false;


	//bool jud_additive_interface_points_inside_Di = false;
	//for (int i = 0; i < temp_V2.rows(); i++) {
	//	for (int j = 0; j < additive_interface_points.size(); j++) {
	//		for (int k = 0; k < additive_interface_points[j].size(); k++) {
	//			//计算additive_interface_points[j][k]与temp_V2[i]的距离
	//			double distance = std::sqrt(pow(additive_interface_points[j][k].x() - temp_V2(i, 0), 2) + pow(additive_interface_points[j][k].y() - temp_V2(i, 1), 2) + pow(additive_interface_points[j][k].z() - temp_V2(i, 2), 2));
	//			if(distance < 0.0000001) {
	//				jud_additive_interface_points_inside_Di = true;
	//				break;
	//			}
	//		}
	//	}
	//}
	//if (jud_additive_interface_points_inside_Di == true)
	//		return false;
	
	is_max_volume = false;
	//Mesh current_mesh;
	sf = clock();
	//CGAL::IO::read_OBJ("MCTS_temp\\" + mesh_target + "\\temp_new_Di.obj", current_mesh);
	ef = clock();
	//sum_file += (double)(ef - sf) / CLOCKS_PER_SEC;
	new_volume = CGAL::Polygon_mesh_processing::volume(new_Di);
	if (new_volume > max_modified_volume) {
		last_max_volume = max_modified_volume;
		max_modified_volume = new_volume;
		modified_cutting_planes_by_MTCS = save_current_modified_cutting_planes;
		is_max_volume = true;
	}
	end_time2 = clock();
	if (flag_final == true) {
		best_modifed_cutting_planes = real_final_cutting_planes;
	}
	//std::cout << "AD time: " << (double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;
	//std::cout << "file time: " << sum_file << "s" << std::endl;
	return true;
}
bool ReFab::Add_triangles_of_cutting_planes(Slicer_2& slicer)
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

void ReFab::Anticlockwise(vector<vector<int>>& real_cutting_plane_triangles, Slicer_2 all_slicer)
{
	for (int t = 0; t < real_cutting_plane_triangles.size(); t++) {
		double max_x = -10000000;
		int index_point;
		for (int i = 0; i < real_cutting_plane_triangles[t].size(); i++)
		{
			if (all_slicer.positions[real_cutting_plane_triangles[t][i]][0] > max_x) {
				max_x = all_slicer.positions[real_cutting_plane_triangles[t][i]][0];
				index_point = i;
			}
		}
		double d = (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point - 1 + real_cutting_plane_triangles[t].size()) % real_cutting_plane_triangles[t].size()]][0])
			* (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point + 1) % real_cutting_plane_triangles[t].size()]][1])
			- (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point - 1 + real_cutting_plane_triangles[t].size()) % real_cutting_plane_triangles[t].size()]][1])
			* (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point + 1) % real_cutting_plane_triangles[t].size()]][0]);
		if (d > 0)
			reverse(real_cutting_plane_triangles[t].begin(), real_cutting_plane_triangles[t].end());
	}
}

//判断两个向量几乎相同
bool ReFab::isEqual(Point3ff v1, Point3ff v2)
{
	if (v1.x - v2.x > -1e-2 && v1.x - v2.x < 1e-2)
		if (v1.y - v2.y > -1e-2 && v1.y - v2.y < 1e-2)
			if (v1.z - v2.z > -1e-2 && v1.z - v2.z < 1e-2)
				return true;
	return false;
}

bool ReFab::isSameDirection(Point3ff v1, Point3ff v2)
{
	if (Dot(v1, v2) > 0)
		return true;
	else
		return false;
}

//在-direction方向下看去,三点是否为逆时针
int ReFab::isAnticlockwise(Point_3 p1, Point_3 p2, Point_3 p3)
{
	double ans = (p2.x() - p1.x()) * (p2.y() - p3.y()) - (p2.y() - p1.y()) * (p2.x() - p3.x());
	//std::cout << ans << '\n';
	if (ans < 0)		//is anticlockwise
		return 1;
	else if (ans == 0)	//is line
		return 2;
	else if (ans > 0)			//is clockwise
		return 3;
}

bool ReFab::checkHaveNoOtherPoint(Point_3 p1, Point_3 p2, Point_3 p3, std::list<Point_3> pointList)
{
	std::list<Point_3>::iterator iter = pointList.begin();
	for (; iter != pointList.end(); iter++)
	{
		Point_3 p = *iter;
		if (p == p1 || p == p2 || p == p3)
			continue;
		//if point p is in the triangle, (p,p1,p2),(p,p2,p3),(p,p3,p1) should be anticlockwise
		if (isAnticlockwise(p, p1, p2) == 1 && isAnticlockwise(p, p2, p3) == 1 && isAnticlockwise(p, p3, p1) == 1)
			return false;
	}
	return true;
}

vector<int> ReFab::poufen(Slicer_2 slicer, vector<int> index_of_points, bool isAnti)
{
	//ni shen zhen:anticlockwise
	//51 and 53 is false
	//bool isAnti = true;
	vector<Point_3> Base;
	for (int i = 0; i < index_of_points.size(); i++) {
		Point_3 temp_point(slicer.positions[index_of_points[i]][0], slicer.positions[index_of_points[i]][1], slicer.positions[index_of_points[i]][2]);
		Base.push_back(temp_point);
	}

	std::list<Point_3> pointList;
	std::list<int> pointIndex;
	std::vector<int> returnIndex;
	int numOfPoint = Base.size();
	if (!isAnti)
	{
		for (int i = numOfPoint - 1; i >= 0; i--)	//can also use list.reverse()
		{
			pointList.push_back(Base[i]);
			pointIndex.push_back(index_of_points[i]);
		}
	}
	else
	{
		for (int i = 0; i < numOfPoint; i++) {
			pointList.push_back(Base[i]);
			pointIndex.push_back(index_of_points[i]);
		}
	}
	//check number of point
	if (pointList.size() < 3)
		return returnIndex;

	std::list<Point_3>::iterator iter = pointList.begin();
	std::list<int>::iterator iter_2 = pointIndex.begin();
	Point_3 p1, p2, p3;
	int t1, t2, t3;
	p1 = *(iter);
	p2 = *(next(iter));
	p3 = *(next(iter, 2));
	t1 = *(iter_2);
	t2 = *(next(iter_2));
	t3 = *(next(iter_2, 2));
	while (pointList.size() >= 3)
	{
		//if(p1,p2,p3) is anticlockwise and no other point in this polygon
		if (isAnticlockwise(p1, p2, p3) == 1 && checkHaveNoOtherPoint(p1, p2, p3, pointList))
		{
			//use face (p1,p2,p3)
			returnIndex.push_back(t1);
			returnIndex.push_back(t2);
			returnIndex.push_back(t3);

			//delete p2
			if (next(iter) == pointList.end()) {
				pointList.erase(pointList.begin());
				pointIndex.erase(pointIndex.begin());
			}

			else {
				pointList.erase(next(iter));
				pointIndex.erase(next(iter_2));
			}

			if (pointList.size() >= 3)
			{
				if (next(iter) == pointList.end())
				{
					p2 = *(pointList.begin());
					p3 = *(++pointList.begin());
					t2 = *(pointIndex.begin());
					t3 = *(++pointIndex.begin());
				}
				else if (next(iter, 2) == pointList.end())
				{
					p2 = *(next(iter));
					p3 = *(pointList.begin());
					t2 = *(next(iter_2));
					t3 = *(pointIndex.begin());
				}
				else
				{
					p2 = *(next(iter));
					p3 = *(next(iter, 2));
					t2 = *(next(iter_2));
					t3 = *(next(iter_2, 2));
				}
			}
		}
		else if (isAnticlockwise(p1, p2, p3) == 2)
		{
			//p1->p2->p3 is a line, so delete p2
			if (next(iter) == pointList.end()) {
				pointList.erase(pointList.begin());
				pointIndex.erase(pointIndex.begin());
			}

			else {
				pointList.erase(next(iter));
				pointIndex.erase(next(iter_2));
			}
			if (pointList.size() >= 3)
			{
				if (next(iter) == pointList.end())
				{
					p2 = *(pointList.begin());
					p3 = *(++pointList.begin());
					t2 = *(pointIndex.begin());
					t3 = *(++pointIndex.begin());
				}
				else if (next(iter, 2) == pointList.end())
				{
					p2 = *(next(iter));
					p3 = *(pointList.begin());
					t2 = *(next(iter_2));
					t3 = *(pointIndex.begin());
				}
				else
				{
					p2 = *(next(iter));
					p3 = *(next(iter, 2));
					t2 = *(next(iter_2));
					t3 = *(next(iter_2, 2));
				}
			}
		}
		else
		{
			iter++;
			iter_2++;
			if (iter == pointList.end())
			{
				iter = pointList.begin();
				p1 = *(iter);
				p2 = *(next(iter));
				p3 = *(next(iter, 2));
				iter_2 = pointIndex.begin();
				t1 = *(iter_2);
				t2 = *(next(iter_2));
				t3 = *(next(iter_2, 2));
			}
			else if (next(iter) == pointList.end())
			{
				p1 = *(iter);
				p2 = *(pointList.begin());
				p3 = *(++pointList.begin());
				t1 = *(iter_2);
				t2 = *(pointIndex.begin());
				t3 = *(++pointIndex.begin());
			}
			else if (next(iter, 2) == pointList.end())
			{
				p1 = *(iter);
				p2 = *(next(iter));
				p3 = *(pointList.begin());
				t1 = *(iter_2);
				t2 = *(next(iter_2));
				t3 = *(pointIndex.begin());
			}
			else
			{
				p1 = *(iter);
				p2 = *(next(iter));
				p3 = *(next(iter, 2));
				t1 = *(iter_2);
				t2 = *(next(iter_2));
				t3 = *(next(iter_2, 2));
			}
		}
	}
	return returnIndex;
}

int ReFab::Calculate_current_num_of_connected_components(double& volume_offset, double& volume_ori, int& total_num_of_complete_planes, vector<Point_3>& barycenter_of_planes, vector<Point_3>& selected_planes_normal, vector<vector<Point_3>>& all_vertices_in_components_input, int interleaving_iteration)
{
	barycenter_of_cutting_planes[1] = barycenter_of_cutting_planes[0];
	barycenter_of_cutting_planes[0] = barycenter_of_cutting_planes[2];
	barycenter_of_cutting_planes.erase(barycenter_of_cutting_planes.begin() + 2);

	int index_last_useful_planes = 0;
	for (int i = 1; i < Selected_Planes.size(); i++)
		if (flag_useful_cutting_planes[i] == true)
			index_last_useful_planes = i; 
	Mesh temp_mesh,temp_mesh2;
	if(interleaving_iteration == 0)
		CGAL::IO::read_OFF("temp_vis\\"+mesh_target+"\\(C)new_D_I-" + to_string(index_last_useful_planes) + ".off", temp_mesh);
	else
		CGAL::IO::read_OBJ("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(Selected_Planes.size()-1) + ".obj", temp_mesh);
	double volume_1 = CGAL::Polygon_mesh_processing::volume(temp_mesh);
	CGAL::IO::read_OBJ("temp_vis\\"+mesh_target+"\\before_offset.obj", temp_mesh2);
	double volume_2 = CGAL::Polygon_mesh_processing::volume(temp_mesh2);
	volume_offset = volume_1;
	volume_ori = volume_2;
	record_data[4] = volume_1 / volume_of_target_model;
	record_data[5] = Selected_Planes.size();
	//barycenter_of_cutting_planes.clear();
	
	//计算每个切平面移除mesh的每个连通分量
	Eigen::MatrixXd temp_V;
	Eigen::MatrixXi temp_F;
	int total_num_of_connected_components = 0;
	cont_components_of_each_planes.resize(Selected_Planes.size());
	index_of_additive_interface.clear();
	vector<int> flag_additive_inter_face_been_coverd(additive_interface_points.size());
	for (int i = 0; i < additive_interface_points.size(); i++)
		flag_additive_inter_face_been_coverd[i] = -1;

	for (int i = 0; i < Selected_Planes.size(); i++) {
		Mesh mesh;
		if (interleaving_iteration == 0) {
			igl::readOBJ("temp_vis\\" + mesh_target + "\\removed_triangles-" + to_string(i) + ".obj", temp_V, temp_F);
			//igl::writeOFF("temp_vis\\" + mesh_target + "\\removed_triangles-" + to_string(i) + ".off", temp_V, temp_F);
			CGAL::IO::read_OBJ("temp_vis\\" + mesh_target + "\\removed_triangles-" + to_string(i) + ".obj", mesh);
		}
		else {
			igl::readOBJ("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(i) + ".obj", temp_V, temp_F);
			//igl::writeOFF("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(i) + ".off", temp_V, temp_F);
			CGAL::IO::read_OBJ("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(i) + ".obj", mesh);
		}
		std::vector<std::vector<int>> connected_components_triangles;
		connected_components_triangles.clear();
		if (i == 0)
			cont_components_of_each_planes[i] = 0;
		else
			cont_components_of_each_planes[i] = cont_components_of_each_planes[i - 1];

		if(flag_useful_cutting_planes[i] == false)
			continue;

		vector<bool> flag_face_been_search(mesh.number_of_faces());
		for (int j = 0; j < mesh.number_of_faces(); j++)
			flag_face_been_search[j] = false;
		typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
		const double bound = std::cos(0.75 * CGAL_PI);
		std::vector<face_descriptor> cc;
		for (face_descriptor fd : faces(mesh))
		{
			cc.clear();
			if (flag_face_been_search[fd.idx()] == true)
				continue;
			PMP::connected_component(fd,
				mesh,
				std::back_inserter(cc));
			std::vector<int> temp_points;
			connected_components_triangles.push_back(temp_points);
			total_num_of_connected_components++;
			cont_components_of_each_planes[i]++;
			selected_planes_normal.push_back(Selected_Planes[i].normal);
			for (int j = 0; j < cc.size(); j++) {
				flag_face_been_search[cc[j].idx()] = true;
				connected_components_triangles[connected_components_triangles.size() - 1].push_back(cc[j].idx());
			}

			//在这里先cut一下，计算其采样点(二维面上固定步长采样)的梯度，然后根据梯度的方向来确定新的切平面。将计算出的新切平面normal的记录下来，用于拓展步骤的启发
			Slicer_2 current_connected_component_mesh;
			current_connected_component_mesh.positions.resize(temp_V.rows());
			current_connected_component_mesh.triangles.resize(connected_components_triangles[connected_components_triangles.size() - 1].size());
			for (int j = 0; j < temp_V.rows(); j++) {
				current_connected_component_mesh.positions[j][0] = temp_V(j, 0);
				current_connected_component_mesh.positions[j][1] = temp_V(j, 1);
				current_connected_component_mesh.positions[j][2] = temp_V(j, 2);
			}
			current_connected_component_mesh.triangles.resize(connected_components_triangles[connected_components_triangles.size() - 1].size());
			for (int j = 0; j < connected_components_triangles[connected_components_triangles.size() - 1].size(); j++) {
				current_connected_component_mesh.triangles[j][0] = temp_F(connected_components_triangles[connected_components_triangles.size() - 1][j], 0);
				current_connected_component_mesh.triangles[j][1] = temp_F(connected_components_triangles[connected_components_triangles.size() - 1][j], 1);
				current_connected_component_mesh.triangles[j][2] = temp_F(connected_components_triangles[connected_components_triangles.size() - 1][j], 2);
			}

			Plane ori_plane;
			ori_plane.normal = Selected_Planes[i].normal;
			//Point_3 current_barycenter = Calculate_Barycenter(current_connected_component_mesh, ori_plane.normal);
			ori_plane.origin = Point_3(barycenter_of_cutting_planes[cont_components_of_each_planes[i] - 1].x(),
				barycenter_of_cutting_planes[cont_components_of_each_planes[i] - 1].y(), barycenter_of_cutting_planes[cont_components_of_each_planes[i] - 1].z());
			//ori_plane.origin = current_barycenter;
			//barycenter_of_cutting_planes.push_back(current_barycenter);

			Slicer_2 temp_connected_component_mesh = current_connected_component_mesh;
			double mean_scalar_value;
			//Eigen::Vector3d expected_normal = Expected_Cutting_Planes_normals_by_gradientField(temp_connected_component_mesh, ori_plane, mean_scalar_value);  //暂时不用梯度场给定期望方向
			Eigen::Vector3d expected_normal;
			expected_normals_of_cutting_planes.push_back(expected_normal);
			mean_scalar_value_of_cutting_planes.push_back(mean_scalar_value);

			//计算当前连通分量包含的增材接触面   
			bool flag_contain_additive_interface;
			vector<int> temp_vec;
			index_of_additive_interface.push_back(temp_vec);

			for (int t = 0; t < coverd_additive_components_by_selected_planes[i].size(); t++) {
				index_of_additive_interface[index_of_additive_interface.size() - 1].push_back(coverd_additive_components_by_selected_planes[i][t]);
			}
			//若某属于同一个selected_cutting_planes的连通分量已覆盖增材接触面，则删除其它连通分量覆盖的接触面
			for (int t = 0; t < coverd_additive_components_by_selected_planes[i].size(); t++) {
				flag_contain_additive_interface = false;
				for (int j = 0; j < connected_components_triangles[connected_components_triangles.size() - 1].size(); j++) {
					for (int k = 0; k < 3; k++) {
						Point_3 temp_point(temp_V(current_connected_component_mesh.triangles[j][k], 0), temp_V(current_connected_component_mesh.triangles[j][k], 1), temp_V(current_connected_component_mesh.triangles[j][k], 2));
						for (int m = 0; m < additive_interface_points[coverd_additive_components_by_selected_planes[i][t]].size(); m++) {
							double distance = std::sqrt(std::pow(temp_point.x() - additive_interface_points[coverd_additive_components_by_selected_planes[i][t]][m].x(), 2) +
								std::pow(temp_point.y() - additive_interface_points[coverd_additive_components_by_selected_planes[i][t]][m].y(), 2) +
								std::pow(temp_point.z() - additive_interface_points[coverd_additive_components_by_selected_planes[i][t]][m].z(), 2));
							if (distance < 0.001) {
								flag_additive_inter_face_been_coverd[coverd_additive_components_by_selected_planes[i][t]] = index_of_additive_interface.size() - 1;
								flag_contain_additive_interface = true;
								break;
							}
						}
						if (flag_contain_additive_interface == true)
							break;
					}
					if (flag_contain_additive_interface == true)
						break;
				}
			}
		
			
			vector<Point_3> temp_all_vertices_in_components;
			for (int j = 0; j < connected_components_triangles[connected_components_triangles.size() - 1].size(); j++) {
				for (int k = 0; k < 3; k++) {
					Point_3 temp_point(temp_V(current_connected_component_mesh.triangles[j][k], 0), temp_V(current_connected_component_mesh.triangles[j][k], 1), temp_V(current_connected_component_mesh.triangles[j][k], 2));
					temp_all_vertices_in_components.push_back(temp_point);
				}
			}
			//Delete duplicate elements in temp_all_vertices_in_components
			for (int j = 0; j < temp_all_vertices_in_components.size(); j++) {
				for (int k = j + 1; k < temp_all_vertices_in_components.size(); k++) {
					if (temp_all_vertices_in_components[j] ==  temp_all_vertices_in_components[k]) {
						temp_all_vertices_in_components.erase(temp_all_vertices_in_components.begin() + k);
						k--;
					}
				}
			}
			all_vertices_in_components_input.push_back(temp_all_vertices_in_components);
		}
	}

	for(int i =0;i<index_of_additive_interface.size();i++)
		for (int j = 0; j < index_of_additive_interface[i].size(); j++) {
			if (flag_additive_inter_face_been_coverd[index_of_additive_interface[i][j]] != i && flag_additive_inter_face_been_coverd[index_of_additive_interface[i][j]] != -1) {
				index_of_additive_interface[i].erase(index_of_additive_interface[i].begin() + j);
				j--;
			}
		}

	//从后往前删除过多考虑的增材接触面
	/*for (int j = index_of_additive_interface.size() - 1; j >= 0; j--) {
		for (int jj = 0; jj < index_of_additive_interface[j].size(); jj++) {
			for (int k = 0; k < j; k++) {
				for (int kk = 0; kk < index_of_additive_interface[k].size(); kk++) {
					if (index_of_additive_interface[j][jj] == index_of_additive_interface[k][kk]) {
						index_of_additive_interface[k].erase(index_of_additive_interface[k].begin() + kk);
						break;
					}
				}
			}
		}
	}*/

	total_num_of_complete_planes = Selected_Planes.size();
	barycenter_of_planes = barycenter_of_cutting_planes;
	all_vertices_in_components = all_vertices_in_components_input;

	Eigen::MatrixXd temp_V2;
	Eigen::MatrixXi temp_F2;
	if (interleaving_iteration == 0) {
		igl::readOBJ("temp_vis\\" + mesh_target + "\\removed_triangles-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::writeOBJ("MCTS_temp\\" + mesh_target + "\\removed_triangles-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::readOBJ("temp_vis\\" + mesh_target + "\\new_D_I-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::writeOBJ("MCTS_temp\\" + mesh_target + "\\new_D_I-" + to_string(0) + ".obj", temp_V2, temp_F2);
	}
	else {
		igl::readOBJ("MCTS_temp\\" + mesh_target + "\\final_removed_triangles-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::writeOBJ("MCTS_temp\\" + mesh_target + "\\removed_triangles-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::readOBJ("MCTS_temp\\" + mesh_target + "\\final_new_D_I-" + to_string(0) + ".obj", temp_V2, temp_F2);
		igl::writeOBJ("MCTS_temp\\" + mesh_target + "\\new_D_I-" + to_string(0) + ".obj", temp_V2, temp_F2);
	}
	return total_num_of_connected_components;
}

double ReFab::distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
	double distance = std::sqrt(std::pow(a.x() - b.x(), 2) +
		std::pow(a.y() - b.y(), 2) +
		std::pow(a.z() - b.z(), 2));
	return distance;
}


double ReFab::trilinearInterpolation(const Eigen::Vector3d& interpolationPoint, const std::vector<DataPoint>& dataPoints)
{
	std::vector<DataPoint> referencePoints(8);
	double min_dis = 1000000;
	Eigen::Vector3d closestPoint;
	double closetValue;

	// 计算离插值点最近的基准点
	for (const auto& dataPoint : dataPoints) {
		double distance = std::sqrt(std::pow(interpolationPoint.x() - dataPoint.point.x(), 2) +
			std::pow(interpolationPoint.y() - dataPoint.point.y(), 2) +
			std::pow(interpolationPoint.z() - dataPoint.point.z(), 2));
		if (distance < min_dis) {
			min_dis = distance;
			closestPoint.x() = dataPoint.point.x();
			closestPoint.y() = dataPoint.point.y();
			closestPoint.z() = dataPoint.point.z();
			closetValue = dataPoint.value;
		}
	}
	if (closestPoint.x() <= interpolationPoint.x() && closestPoint.y() <= interpolationPoint.y() && closestPoint.z() <= interpolationPoint.z())  //1
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z() + minimal_d);
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() + minimal_d, closestPoint.z() + minimal_d);
	}
	else if (closestPoint.x() >= interpolationPoint.x() && closestPoint.y() <= interpolationPoint.y() && closestPoint.z() <= interpolationPoint.z())  //2
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() + minimal_d, closestPoint.z() + minimal_d);
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z() + minimal_d);
	}
	else if (closestPoint.x() <= interpolationPoint.x() && closestPoint.y() >= interpolationPoint.y() && closestPoint.z() <= interpolationPoint.z()) //3
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z() + minimal_d);
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() - minimal_d, closestPoint.z() + minimal_d);
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z() + minimal_d);
	}
	else if (closestPoint.x() >= interpolationPoint.x() && closestPoint.y() >= interpolationPoint.y() && closestPoint.z() <= interpolationPoint.z())  //4
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() - minimal_d, closestPoint.z() + minimal_d);
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z() + minimal_d);
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z() + minimal_d);
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() + minimal_d);
	}
	else if (closestPoint.x() <= interpolationPoint.x() && closestPoint.y() <= interpolationPoint.y() && closestPoint.z() >= interpolationPoint.z())  //5
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z() - minimal_d);
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() + minimal_d, closestPoint.z() - minimal_d);
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() + minimal_d, closestPoint.z());
	}
	else if (closestPoint.x() >= interpolationPoint.x() && closestPoint.y() <= interpolationPoint.y() && closestPoint.z() >= interpolationPoint.z())  //6
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() + minimal_d, closestPoint.z() - minimal_d);
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z() - minimal_d);
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() + minimal_d, closestPoint.z());
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() + minimal_d, closestPoint.z());
	}
	else if (closestPoint.x() <= interpolationPoint.x() && closestPoint.y() >= interpolationPoint.y() && closestPoint.z() >= interpolationPoint.z())  //7
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z() - minimal_d);
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() - minimal_d, closestPoint.z() - minimal_d);
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x() + minimal_d, closestPoint.y(), closestPoint.z());
	}
	else if (closestPoint.x() >= interpolationPoint.x() && closestPoint.y() >= interpolationPoint.y() && closestPoint.z() >= interpolationPoint.z())  //8
	{
		referencePoints[0].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() - minimal_d, closestPoint.z() - minimal_d);
		referencePoints[1].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z() - minimal_d);
		referencePoints[2].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[3].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z() - minimal_d);
		referencePoints[4].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[5].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y() - minimal_d, closestPoint.z());
		referencePoints[6].point = Eigen::Vector3d(closestPoint.x() - minimal_d, closestPoint.y(), closestPoint.z());
		referencePoints[7].point = Eigen::Vector3d(closestPoint.x(), closestPoint.y(), closestPoint.z());
	}
	
	for (int i = 0; i < 8; i++) {
		bool flag_find_point = false;
		for (const auto& dataPoint : dataPoints) {   //这里需要加速
			if (distance(dataPoint.point, referencePoints[i].point) <= 0.0001) {
				referencePoints[i].value = dataPoint.value;
				flag_find_point = true;
				break;
			}	
		}
		if (flag_find_point == false) {   //边界不连续
			referencePoints[i].value = closetValue;
		}	
	}

	/*for (int i = 0; i < 8; i++)
		cout << i << " " << referencePoints[i].point.x() << " " << referencePoints[i].point.y() << " " << referencePoints[i].point.z() << " " << referencePoints[i].value << endl;
	cout << endl;*/

	//三线性插值
	double x = interpolationPoint.x();
	double y = interpolationPoint.y();
	double z = interpolationPoint.z();
	double label000 = referencePoints[0].value;
	double label001 = referencePoints[4].value;
	double label010 = referencePoints[2].value;
	double label011 = referencePoints[6].value;
	double label100 = referencePoints[1].value;
	double label101 = referencePoints[5].value;
	double label110 = referencePoints[3].value;
	double label111 = referencePoints[7].value;

	double xd = (x - referencePoints[0].point.x()) / (referencePoints[1].point.x() - referencePoints[0].point.x());
	double yd = (y - referencePoints[0].point.y()) / (referencePoints[2].point.y() - referencePoints[0].point.y());
	double zd = (z - referencePoints[0].point.z()) / (referencePoints[4].point.z() - referencePoints[0].point.z());

	// Performing trilinear interpolation
	double c00 = label000 * (1 - xd) + label100 * xd;
	double c01 = label001 * (1 - xd) + label101 * xd;
	double c10 = label010 * (1 - xd) + label110 * xd;
	double c11 = label011 * (1 - xd) + label111 * xd;

	double c0 = c00 * (1 - yd) + c10 * yd;
	double c1 = c01 * (1 - yd) + c11 * yd;

	double c = c0 * (1 - zd) + c1 * zd;

	return c;
	
}

Eigen::Vector3d ReFab::gradientField(const std::vector<DataPoint>& vertices, const Eigen::Vector3d& p, double delta) 
{
	double x = p.x();
	double y = p.y();
	double z = p.z();

	double label_xp = trilinearInterpolation({ x + delta, y, z },vertices);
	double label_xn = trilinearInterpolation({ x - delta, y, z },vertices);
	double label_yp = trilinearInterpolation({ x, y + delta, z },vertices);
	double label_yn = trilinearInterpolation({ x, y - delta, z },vertices);
	double label_zp = trilinearInterpolation({ x, y, z + delta },vertices);
	double label_zn = trilinearInterpolation({ x, y, z - delta },vertices);

	double grad_x = (label_xp - label_xn) / (2 * delta);
	double grad_y = (label_yp - label_yn) / (2 * delta);
	double grad_z = (label_zp - label_zn) / (2 * delta);

	Eigen::Vector3d gradient(grad_x, grad_y, grad_z);
	gradient.normalize();
	return gradient;
}

Eigen::Vector3d ReFab::Expected_Cutting_Planes_normals_by_gradientField(Slicer_2 slicer, Plane plane, double& mean_scalar_value)
{
	typedef CGAL::Simple_cartesian<double> Kernel;
	typedef Kernel::Point_2 Point_2;
	typedef CGAL::Polygon_2<Kernel> Polygon_2;

	int num_triangles_before_cutting = slicer.positions.size();
	slicer.normal[0] = plane.normal.x();
	slicer.normal[1] = plane.normal.y();
	slicer.normal[2] = plane.normal.z();
	slicer.origin[0] = plane.origin.x();
	slicer.origin[1] = plane.origin.y();
	slicer.origin[2] = plane.origin.z();
	slicer.cut();
	if (num_triangles_before_cutting == slicer.positions.size())
		cout << "no cutting" << endl;

	vector<Vector3> all_the_points_of_cutting_planes;
	std::map<std::pair<int, int>, int>::iterator it;
	for (it = slicer.intersections.begin(); it != slicer.intersections.end(); it++)
		all_the_points_of_cutting_planes.push_back(Vector3(slicer.positions[it->second][0], slicer.positions[it->second][1], slicer.positions[it->second][2]));
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(plane.normal.x(), plane.normal.y(), plane.normal.z());
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	vector<Vector3> temp_all_the_points_of_cutting_planes = all_the_points_of_cutting_planes;
	for (int j = 0; j < temp_all_the_points_of_cutting_planes.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = temp_all_the_points_of_cutting_planes[j].m_x;
		temp_V(1, 0) = temp_all_the_points_of_cutting_planes[j].m_y;
		temp_V(2, 0) = temp_all_the_points_of_cutting_planes[j].m_z;
		temp_V = rotMatrix.inverse() * temp_V;
		temp_all_the_points_of_cutting_planes[j] = Vector3(temp_V(0, 0), temp_V(1, 0), temp_V(2, 0));
	}
	
	Vector3 centroid = { 0, 0 ,0};
	for (const  Vector3& p : temp_all_the_points_of_cutting_planes) {
		centroid.m_x += p.m_x;
		centroid.m_y += p.m_y;
	}
	centroid.m_x /= temp_all_the_points_of_cutting_planes.size();
	centroid.m_y /= temp_all_the_points_of_cutting_planes.size();
	int n = temp_all_the_points_of_cutting_planes.size();
	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < n - i - 1; ++j) {
			double angle1 = atan2(temp_all_the_points_of_cutting_planes[j].m_y - centroid.m_y, temp_all_the_points_of_cutting_planes[j].m_x - centroid.m_x);
			double angle2 = atan2(temp_all_the_points_of_cutting_planes[j+1].m_y - centroid.m_y, temp_all_the_points_of_cutting_planes[j+1].m_x - centroid.m_x);
			if (angle1 > angle2) {
				std::swap(temp_all_the_points_of_cutting_planes[j], temp_all_the_points_of_cutting_planes[j + 1]);
			}
		}
	}
	vector<Point_2> points_2d;
	for(int i =0;i<temp_all_the_points_of_cutting_planes.size();i++)
		points_2d.push_back(Point_2(temp_all_the_points_of_cutting_planes[i].m_x, temp_all_the_points_of_cutting_planes[i].m_y));
	Polygon_2 polygon(points_2d.begin(), points_2d.end());

	//在多边形内采样
	std::vector<Point_2> sampling_points_in_cutting_planes;
	CGAL::Bbox_2 bbox = polygon.bbox();
	Point_2 out_minC(bbox.xmin(), bbox.ymin()), out_maxC(bbox.xmax(), bbox.ymax());
	double x(out_minC[0]);
	float sampling_step = 1.0;
	while (x < out_maxC[0])
	{
		double y(out_minC[1]);
		while (y < out_maxC[1])
		{
			Point_2 p(x, y);
			if (polygon.bounded_side(p) == CGAL::ON_BOUNDED_SIDE)
			{
				sampling_points_in_cutting_planes.push_back(p);
			}
			y += sampling_step;
		}
		x += sampling_step;
	}

	//旋转回去
	rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vectorBefore).toRotationMatrix();
	vector<Eigen::Vector3d> temp_sampling_points_in_cutting_planes;
	for (int i = 0; i < sampling_points_in_cutting_planes.size(); i++)
		temp_sampling_points_in_cutting_planes.push_back(Eigen::Vector3d(sampling_points_in_cutting_planes[i].x(),sampling_points_in_cutting_planes[i].y(), temp_all_the_points_of_cutting_planes[0].m_z));
	for (int j = 0; j < temp_sampling_points_in_cutting_planes.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = temp_sampling_points_in_cutting_planes[j].x();
		temp_V(1, 0) = temp_sampling_points_in_cutting_planes[j].y();
		temp_V(2, 0) = temp_sampling_points_in_cutting_planes[j].z();
		temp_V = rotMatrix.inverse() * temp_V;
		temp_sampling_points_in_cutting_planes[j] = Eigen::Vector3d(temp_V(0, 0), temp_V(1, 0), temp_V(2, 0));
	}
	//ofstream off("aaaaa.obj");
	////std::map<std::pair<int, int>, int>::iterator it;
	//for (int i =0;i< temp_sampling_points_in_cutting_planes.size();i++)
	//		off << "v " << temp_sampling_points_in_cutting_planes[i].x()<< " " << temp_sampling_points_in_cutting_planes[i].y() << " " << temp_sampling_points_in_cutting_planes[i].z()<< endl;

	//计算所有采样点在梯度场中的梯度
	Eigen::Vector3d new_normal(0,0,0);
	double current_mean_scalar_value = 0;
	//ofstream test_gradient("gradient.txt");
	//ofstream test_interpolation("test_interpolation.txt");
	for(double i=0; i< temp_sampling_points_in_cutting_planes.size();i++)
	{
		double interpolatedValue = trilinearInterpolation(temp_sampling_points_in_cutting_planes[i], current_green_points);
		Eigen::Vector3d gradient = gradientField(current_green_points, temp_sampling_points_in_cutting_planes[i], 0.001);
		//test_interpolation << temp_sampling_points_in_cutting_planes[i].x() << " " << temp_sampling_points_in_cutting_planes[i].y() << " " << temp_sampling_points_in_cutting_planes[i].z() << " " << interpolatedValue << endl;
		//test_gradient << temp_sampling_points_in_cutting_planes[i].x() << " " << temp_sampling_points_in_cutting_planes[i].y() << " " << temp_sampling_points_in_cutting_planes[i].z() << " " << gradient.x() << " " << gradient.y() << " " << gradient.z() << endl;
		current_mean_scalar_value += interpolatedValue;
		new_normal += gradient;
	}
	current_mean_scalar_value /= temp_sampling_points_in_cutting_planes.size();
	mean_scalar_value = current_mean_scalar_value;
	new_normal /= temp_sampling_points_in_cutting_planes.size();
	new_normal.normalize();

	///////////////用于测试插值，之后注释///////////////
	///*for(float x = 11;x< 23;x+=1)
	//	for(float y =-2;y< 10;y+=1)
	//		for(float z=75;z<87;z+=1)
	//			test_sampling_points.push_back(Eigen::Vector3d(x, y, z));*/

	//ofstream test_gradient("gradient.txt");
	//ofstream test_interpolation("test_interpolation.txt");
	////ofstream test_obj("test.obj");
	//for (double i = 0; i < test_sampling_points.size(); i++)
	//{
	//	double interpolatedValue = trilinearInterpolation(test_sampling_points[i], current_green_points);
	//	Eigen::Vector3d gradient = gradientField(current_green_points, test_sampling_points[i], 0.001);
	//	test_interpolation << test_sampling_points[i].x() << " " << test_sampling_points[i].y() << " " << test_sampling_points[i].z() << " " << interpolatedValue << endl;
	//	test_gradient << test_sampling_points[i].x() << " " << test_sampling_points[i].y() << " " << test_sampling_points[i].z() << " " << gradient.x() << " " << gradient.y() << " " << gradient.z() << endl;
	//	//test_obj << "v " << test_sampling_points[i].x() << " " << test_sampling_points[i].y() << " " << test_sampling_points[i].z() << endl;
	//}
	//////////////////////////////////////////

	//ofstream off("aaaaa.obj");
	////std::map<std::pair<int, int>, int>::iterator it;
	//for (int i =0;i< temp_all_the_points_of_cutting_planes.size();i++)
	//		off << "v " << temp_all_the_points_of_cutting_planes[i].m_x << " " << temp_all_the_points_of_cutting_planes[i].m_y << " "  << temp_all_the_points_of_cutting_planes[i].m_z  << endl;

	return new_normal;
}

Point_3 ReFab::Calculate_Barycenter(Slicer_2 current_connected_component_mesh, Point_3 normal)
{
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Vector3d vector_after(normal.x(), normal.y(), normal.z());
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int j = 0; j < current_connected_component_mesh.positions.size(); j++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = current_connected_component_mesh.positions[j][0];
		temp_V(1, 0) = current_connected_component_mesh.positions[j][1];
		temp_V(2, 0) = current_connected_component_mesh.positions[j][2];
		temp_V = rotMatrix.inverse() * temp_V;
		current_connected_component_mesh.positions[j][0] = temp_V(0, 0);
		current_connected_component_mesh.positions[j][1] = temp_V(1, 0);
		current_connected_component_mesh.positions[j][2] = temp_V(2, 0);
	}

	//计算z值最小的current_connected_component_mesh.position
	double min_z = 100000;
	for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if(current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][2] < min_z)
				min_z = current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][2];
		}
	}
		
	ofstream file("ttttt.obj");
	for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			file << "v " << current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][0] << " "
				<< current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][1] <<" "
				<<current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][2]<<endl;
		}
	}
	Slicer_2 selected_triangles_from_slicer;
	selected_triangles_from_slicer.clear();
	for (int j = 0; j < current_connected_component_mesh.triangles.size(); j++) {
		for (int k = 0; k < 3; k++) {
			if (abs(current_connected_component_mesh.positions[current_connected_component_mesh.triangles[j][k]][2] - min_z) < 0.01) {
				selected_triangles_from_slicer.triangles.push_back(current_connected_component_mesh.triangles[j]);
				current_connected_component_mesh.triangles.erase(current_connected_component_mesh.triangles.begin() + j);
				j--;
				break;
			}
		}
	}
	for (int j = 0; j < selected_triangles_from_slicer.triangles.size(); j++) {
		selected_triangles_from_slicer.positions.push_back(current_connected_component_mesh.positions[selected_triangles_from_slicer.triangles[j][0]]);
		selected_triangles_from_slicer.positions.push_back(current_connected_component_mesh.positions[selected_triangles_from_slicer.triangles[j][1]]);
		selected_triangles_from_slicer.positions.push_back(current_connected_component_mesh.positions[selected_triangles_from_slicer.triangles[j][2]]);
	}

	//Delete duplicate elements in temp_all_vertices_in_components
	for (int j = 0; j < selected_triangles_from_slicer.positions.size(); j++) {
		for (int k = j + 1; k < selected_triangles_from_slicer.positions.size(); k++) {
			if (selected_triangles_from_slicer.positions[j] == selected_triangles_from_slicer.positions[k]) {
				selected_triangles_from_slicer.positions.erase(selected_triangles_from_slicer.positions.begin() + k);
				k--;
			}
		}
	}

	//calculate the barycenter
	Eigen::MatrixXd temp_barycenter;
	temp_barycenter.resize(3, 1);
	double barycenter_x = 0, barycenter_y = 0, Area = 0;
	for (int j = 0; j < selected_triangles_from_slicer.positions.size()-1; j++) {
		double x_1 =selected_triangles_from_slicer.positions[j][0];
		double y_1 =selected_triangles_from_slicer.positions[j][1];
		double x_2 =selected_triangles_from_slicer.positions[j+1][0];
		double y_2 =selected_triangles_from_slicer.positions[j+1][1];
		barycenter_x += (x_1 + x_2) * (x_1 * y_2 - x_2 * y_1);
		barycenter_y += (y_1 + y_2) * (x_1 * y_2 - x_2 * y_1);
		Area += x_1 * y_2 - x_2 * y_1;
	}
	Area /= 2;
	barycenter_x /= 6 * Area;
	barycenter_y /= 6 * Area;
	temp_barycenter(0, 0) = barycenter_x;
	temp_barycenter(1, 0) = barycenter_y;
	temp_barycenter(2, 0) = min_z +0.1;  
	Eigen::Matrix3d rotMatrix_2 = Eigen::Quaterniond::FromTwoVectors(vector_after, vectorBefore).toRotationMatrix();
	temp_barycenter = rotMatrix_2.inverse() * temp_barycenter;
	Point_3 current_barycenter(temp_barycenter(0, 0), temp_barycenter(1, 0), temp_barycenter(2, 0));
	//barycenter_of_cutting_planes.push_back(Point_3(temp_barycenter(0,0), temp_barycenter(1, 0), temp_barycenter(2, 0)));

	return current_barycenter;
}