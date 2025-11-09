
#pragma once

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include"Vector3.h"
#include"visualization.h"
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include<igl/readOFF.h>
#include "MeshPlaneIntersect.hpp"
#include"earcut.hpp"
#include "slicer.h"
//#include"IState.h"
#include <vector>
#include <string>
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
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include "Point3f.h"
#include<direct.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include"Element.h"
#include <unordered_set>
#include<random>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/opengl/glfw/Viewer.h>
#include <glpk.h>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef MeshPlaneIntersect<double, int> Intersector;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Property_map<vertex_descriptor, EK::Point_3> Exact_point_map;
typedef K::Point_3     Point_3;



namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;
#define PI 3.1415926

struct inaccessible_area
{
    int id_to_point;
    int id_ori;
    inaccessible_area(int a = 0, int b = 0) {
        id_to_point = a;
        id_ori = b;
    }
};

struct nozzle
{
    double upper_surface_r;
    double lowwer_surface_r;
    double nozzle__H_total;
    double nozzle_H_half;
};

struct cutter
{
    double cylinder_height;
    double ball_r;
    double carriage_r;
    double carriage_height;
};

struct Exact_vertex_point_map
{
    // typedef for the property map
    typedef boost::property_traits<Exact_point_map>::value_type value_type;
    typedef boost::property_traits<Exact_point_map>::reference reference;
    typedef boost::property_traits<Exact_point_map>::key_type key_type;
    typedef boost::read_write_property_map_tag category;
    // exterior references
    Exact_point_map exact_point_map;
    Mesh* tm_ptr;
    // Converters
    CGAL::Cartesian_converter<K, EK> to_exact;
    CGAL::Cartesian_converter<EK, K> to_input;
    Exact_vertex_point_map()
        : tm_ptr(nullptr)
    {}
    Exact_vertex_point_map(const Exact_point_map& ep, Mesh& tm)
        : exact_point_map(ep), tm_ptr(&tm)
    {
        for (Mesh::Vertex_index v : vertices(tm))
            exact_point_map[v] = to_exact(tm.point(v));
    }
    friend
        reference get(const Exact_vertex_point_map& map, key_type k)
    {
        CGAL_precondition(map.tm_ptr != nullptr);
        return map.exact_point_map[k];
    }
    friend
        void put(const Exact_vertex_point_map& map, key_type k, const EK::Point_3& p)
    {
        CGAL_precondition(map.tm_ptr != nullptr);
        map.exact_point_map[k] = p;
        // create the input point from the exact one
        map.tm_ptr->point(k) = map.to_input(p);
    }
};

struct Vector_pmap_wrapper
{
    std::vector<bool>& vect;
    Vector_pmap_wrapper(std::vector<bool>& v) : vect(v) {}
    friend bool get(const Vector_pmap_wrapper& m, face_descriptor f)
    {
        return m.vect[f];
    }
    friend void put(const Vector_pmap_wrapper& m, face_descriptor f, bool b)
    {
        m.vect[f] = b;
    }
};

struct Plane
{
    Point_3 origin;
    Point_3 normal;
    Plane(Point_3 nor, Point_3 ori) {
        normal = nor;
        origin = ori;
	}
    Plane() {};
};

struct DataPoint {
    Eigen::Vector3d point;
    double value;
    DataPoint() { ; }
    DataPoint(double a, double b, double c, double d) {
        Eigen::Vector3d temp(a, b, c);
        point = temp;
		value = d;
    }
};

struct Point2d {
    double x, y;
    Point2d(double a, double b) {
        x = a;
        y = b;
    }
};

class ReFab {
    public:
    //根目录地址
    std::string rootPath;
    ReFab() {}

    //构造函数
    ReFab(std::string& rootpath, cutter Cutter, nozzle Nozzle)
        :rootPath(rootpath)
    {
        this->Cutter = Cutter;
		this->Nozzle = Nozzle;
    }
    
    ~ReFab() {}
    //计算两个mesh的差
    void Difference_two_mesh(std::string& input1, std::string& input2, std::string output = "differece.off", int type = 0);
    //计算两个mesh的交
    void Intersection_two_mesh(std::string &input1, std::string& input2, std::string output = "intersection.off",  bool type=0);
    //计算Di的基底
    void DetermineBase(string path);
    //计算两个mesh的并
    void Union_two_mesh(std::string& input1, std::string& input2, std::string output = "union.off");
    //对模型进行求交，并进行分割，返回分割后的模型
    void Intersection_and_Split(std::string& input1, std::string& input2, std::string output1 = "output\\interface_s.off",
        std::string output2 = "output\\interface_a.off", std::string output = "output\\intersection.off",bool type=0);
    //对于给定mesh，在其内部均匀采样
    void Mesh_Regular_Sampling_C1(std::string& outside_path, const double& d);
    void Get_sampling_points_in_DS_and_DA(std::string DS_file, std::string DA_file);
    void OrientationSamplePoints();
    std::vector<Point_3> OrientationSamplePoints_normal();
    std::vector<Point_3> OrientationSamplePoints_normal_2();
    void Get_additive_inter_face(Mesh inter, Mesh mesh_a, Mesh mesh_s, std::string output_a, double& sum_additive_interface_area, double& min_z_of_additive_interface);
    double Calculate_area_need_support(Eigen::MatrixXd V_1, Eigen::MatrixXi F_1, double& sum_area);
    void Get_additive_inter_face(std::string input_inter, std::string input_a, std::string input_s, std::string output_a, double& sum_additive_interface_area, double& min_z_of_additive_interface);
    void ReOrientation(std::string& meshname1, std::string& meshname2, bool type);
    void ReOrientation_combined_optimization(std::string& meshname1, std::string& meshname2, bool type);
    Point_3 calculate_normal(Point_3 a, Point_3 b, Point_3 c);
    void Collision_Detection_For_DS_and_DA();
    void Update_scalar_field_DS_DA(std::vector<int> index_removed_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration);
    void Update_scalar_field_DS_DA_with_removed_points(std::vector<int> index_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration);
    void Update_scalar_field_DS_DA_2(std::vector<int> index_removed_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA,int iteration);
    void Update_scalar_field_DS_DA_with_removed_points_2(std::vector<int> index_Di_points, std::vector<Eigen::MatrixXd>& vis_red_points_DS, std::vector<Eigen::MatrixXd>& vis_red_points_DA, int iteration);
    void Greedy_one_by_one_delete_points_of_Di();
    double Cut_mesh_by_one_plane(Plane current_plane, double sum_value_of_green);
    vector<Plane> Sort_cutting_planes(vector<Plane> Selected_Planes, std::vector<std::vector<Point_3>>& additive_interface_points, std::vector<std::vector<int>>& coverd_additive_components_by_selected_planes);
    void Initialize_Cutting_Planes(string additive_interface);
    bool Check_is_one_component(Mesh new_Di);
    void Calculate_mean_curvature(string additive_interface, bool vis, double& mean_mean_principal_curvature, double& abs_mean_mean_principal_curvature);
    void Modify_Cutting_Planes();
    bool Modify_Cutting_Planes_2(vector<Plane> parent_planes, double& new_volume, vector<Plane>& best_modifed_cutting_planes, double& max_modified_volume, double& last_max_volume, vector<Point_3>& new_current_removed_points, bool flag_final, bool& is_max_volume);
    bool Add_triangles_of_cutting_planes(Slicer_2& slicer);
    double Calculate_distance_plane_and_point(Point_3 normal, Point_3 origin, Point_3 p);
    int  Calculate_current_num_of_connected_components(double& volume_diff, double& volume_ori, int& total_num_of_complete_planes, vector<Point_3>& barycenter_of_planes, vector<Point_3>& selected_planes_normal, vector<vector<Point_3>>& all_vertices_in_components, int interleaving_iteration);

    void Anticlockwise(vector<vector<int>>& real_cutting_plane_triangles, Slicer_2 all_slicer);
    vector<int> poufen(Slicer_2 slicer, vector<int> index_of_points, bool isAnti);
    bool checkHaveNoOtherPoint(Point_3 p1, Point_3 p2, Point_3 p3, std::list<Point_3> pointList);
    bool isSameDirection(Point3ff v1, Point3ff v2);
    bool isEqual(Point3ff v1, Point3ff v2);
    int isAnticlockwise(Point_3 p1, Point_3 p2, Point_3 p3);
    double distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b);
    double trilinearInterpolation(const Eigen::Vector3d& interpolationPoint, const std::vector<DataPoint>& dataPoints);
    Eigen::Vector3d gradientField(const std::vector<DataPoint>& vertices, const Eigen::Vector3d& p, double delta);
    Eigen::Vector3d Expected_Cutting_Planes_normals_by_gradientField(Slicer_2 slicer, Plane plane, double& mean_scalar_value);

    void CalculateHangInArea(string input_DI);
    void Modify_DS_and_DA(string input_DS, string input_DA, string input_hang_in_area);
    Point_3 Calculate_Barycenter(Slicer_2 current_connected_component_mesh, Point_3 normal);

    void solveWeightedSetCover(glp_prob* lp,int n,int m,std::vector<int> weights, std::vector<std::set<int>> subsets);

    cutter Cutter;
    nozzle Nozzle;
    std::vector<int> index_base_of_Di;
    std::vector<Vec3> vertices_base;
    std::vector<Point_3> normal_D_S;
    std::vector<Point_3> normal_D_A;
    std::vector<Point_3> sampling_points_in_D_S;
    std::vector<Point_3> sampling_points_in_D_A;   
    std::vector<Point_3> sampling_points_in_D_I;
    std::vector<Point_3> sample_points_sphere;
    std::vector<Eigen::Vector3d> test_sampling_points;
    std::vector<std::vector<inaccessible_area>> all_the_inaccessible_area_DS, all_the_inaccessible_area_DA;
    std::vector<std::vector<inaccessible_area>> all_the_covering_points_and_S, all_the_covering_points_and_A;
    std::vector<std::vector<inaccessible_area>> all_the_covering_points_and_S_2, all_the_covering_points_and_A_2;
    std::map<int, int> map_index_and_DS_points;
    std::map<int, int> map_index_and_DS_points_inv;
    std::map<int, int> map_covering_points_and_DS_points;
    std::map<int, int> map_covering_points_and_DS_points_inv;
    std::map<int, int> map_index_and_DA_points;
    std::map<int, int> map_index_and_DA_points_inv;
    std::map<int, int> map_covering_points_and_DA_points;
    std::map<int, int> map_covering_points_and_DA_points_inv;
    std::vector<std::vector<Point_3>> additive_interface_points;
    std::vector<bool> flag_accessible_points_S, flag_accessible_points_A;
    std::vector<double> color_map;
    vector<DataPoint> current_green_points;
    
    std::vector<bool> flag_use_selected_planes;
    std::vector<Plane> Selected_Planes;
    std::vector<std::vector<int>> coverd_additive_components_by_selected_planes;
    vector<int> index_triangles_of_cutting_planes;
    vector<Point_3> barycenter_of_cutting_planes;
    vector<int> cont_components_of_each_planes;
    /*std::vector<Plane> Current_Planes;
    std::vector<std::vector<Point_3>> Polygons_of_Current_Planes;*/
    
    bool open_vis_red_points;
    bool open_vis_green_points;
    bool* flag_useful_cutting_planes;
    vector<Plane> modified_cutting_planes_by_MTCS;
    string mesh_target;
    double min_z_of_DI;
    double minimal_d;
    vector<Eigen::Vector3d> expected_normals_of_cutting_planes;
    vector<double> mean_scalar_value_of_cutting_planes;
    double volume_of_Di;
    double volume_of_target_model;
    vector<vector<int>> index_of_additive_interface;
    vector<vector<Point_3>> all_vertices_in_components;

    vector<Plane>final_cutting_planes;
    vector<Point_3> points_in_cutting_planes;

    Slicer_2 saved_removed_triangles, saved_new_D_I;
    vector<float> record_data;
};

//找到给定三角形网格（Polyhedron类型）中所有顶点的最大坐标值
double max_coordinate(const Mesh& mesh);

// 辅助函数：判断面是否属于指定的网格
// 通过face的顶点坐标和mesh的位置关系判断
bool is_face_in_mesh(const Mesh& mesh, const Mesh::Face_index& face);






