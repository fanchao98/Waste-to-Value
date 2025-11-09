#pragma once
//Í·ÎÄ¼þ
#include"ReFab.h"
#include<igl\writeSTL.h>
#include"katana.h"
#include<cstdio>

struct My_node
{
    int indegree;
    vector<pair<int, int>> edges;
    My_node()
    {
        indegree = 0;
        edges.clear();
    }
};

class PPS {
public:
    void Calculate_additive_blocks(string meshname2, vector<Plane>final_cutting_planes, vector<vector<Point_3>>points_in_cutting_planes);
    double Calculate_area_need_support_2(Slicer_2& slicer, double& sum_area, bool& jud_flat_area);
    void Segmentation_for_each_removed_mesh(vector<Plane> final_cutting_planes);
    bool Add_triangles_of_cutting_planes(Slicer_2& slicer);
    bool Add_triangles_of_cutting_planes_2(Slicer_2& slicer);
    bool Check_cutting_bottom(Slicer_2& slicer);
    void OrientationSamplePoints();
    double Calculate_area_need_support(Slicer_2& slicer, double& sum_area);
    void Sort_slice_layers();

    string mesh_target;
    std::vector<Point_3> sample_points_sphere;

    vector<vector<vector<Vertex>>> all_slice_points;
    vector<vector<Vertex>> all_slice_normal;
    vector<pair<int, int>> all_index_slice;
    vector<pair<int, int>> new_slicer_order;

    nozzle Nozzle;
};