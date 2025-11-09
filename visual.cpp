#include "visual.h"
#include <direct.h>


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

void Visual::generateModelForRendering(std::vector<std::vector<std::vector<Eigen::MatrixXd>>> lines, string file_name)
{
    float radius = 0.7;
    General_Mesh mesh1;
    mesh1.r = rand() / double(RAND_MAX);
    mesh1.g = rand() / double(RAND_MAX);
    mesh1.b = rand() / double(RAND_MAX);
    std::vector<Vec3> points;
    for (int i = 0; i < lines.size(); i++)
    {
        for (int j = 0; j < lines[i].size(); j++)
        {
            for (int k = 0; k < lines[i][j].size(); k++) {
                Vec3 the_point;
                the_point.m_x = lines[i][j][k](0,0);
                the_point.m_y = lines[i][j][k](1,0);
                the_point.m_z = lines[i][j][k](2,0);
                points.push_back(the_point);
            }
            insert_Line(mesh1, points, radius);
            points.clear();
        } 
    }
    mesh1.genResultMesh((".\\Results\\" + file_name + ".obj").c_str());
    points.clear();
}

void Visual::generateModelForRendering_2(Eigen::Vector3d lines, string file_name)
{
    float radius = 0.1;
    General_Mesh mesh1;
    mesh1.r = rand() / double(RAND_MAX);
    mesh1.g = rand() / double(RAND_MAX);
    mesh1.b = rand() / double(RAND_MAX);
    std::vector<Vec3> points;

    Vec3 the_point;
    the_point.m_x = lines.x();
    the_point.m_y = lines.y();
    the_point.m_z = lines.z();
    Vec3 the_point_2(0, 0, 0);
    points.push_back(the_point_2);
    points.push_back(the_point);

    insert_Line(mesh1, points, radius);
    points.clear();


    mesh1.genResultMesh((".\\Results\\" + file_name + ".obj").c_str());
    points.clear();
}


void Visual::generateModelForRendering_3(Eigen::Vector3d lines, string file_name, std::vector<Eigen::MatrixXd> vis_points)
{
    float radius = 3;
    General_Mesh mesh1;
    mesh1.r = 0.7;
    mesh1.g = 0.15;
    mesh1.b = 0.15;
    std::vector<Vec3> points;

    Vec3 the_point;
    
    Vec3 the_point_2;
    for (int i = 0; i < 1; i++) {
        the_point_2.m_x = vis_points[i](0, 0);
        the_point_2.m_y = vis_points[i](1, 0);
        the_point_2.m_z = vis_points[i](2, 0);

        the_point.m_x = the_point_2.m_x + lines.x() * 100;
        the_point.m_y = the_point_2.m_y + lines.y() * 100;
        the_point.m_z = the_point_2.m_z + lines.z() * 100;
        points.push_back(the_point_2);
        points.push_back(the_point);
        insert_Line(mesh1, points, radius);
        points.clear();
    }
    mesh1.genResultMesh((".\\Results\\" + file_name + ".obj").c_str());
    points.clear();
}

void Visual::generateModelForRendering_4(std::vector<std::vector<Eigen::MatrixXd>> lines, string file_name, vector<vector<int>> all_blocks)
{
    float radius = 0.7;
    std::vector<Vec3> points;
    string derectory = ".\\Results\\" + file_name;
    _mkdir(derectory.c_str());
    for (int block = 0; block < all_blocks.size(); block++) {
        General_Mesh mesh1;
        mesh1.r = rand() / double(RAND_MAX);
        mesh1.g = rand() / double(RAND_MAX);
        mesh1.b = rand() / double(RAND_MAX);
        string file_block = "block" + to_string(block);
        for (int i = 0; i < all_blocks[block].size(); i++) {
            for (int j = 0; j < lines[all_blocks[block][i]].size(); j++) {
                Vec3 the_point;
                the_point.m_x = lines[all_blocks[block][i]][j](0, 0);
                the_point.m_y = lines[all_blocks[block][i]][j](1, 0);
                the_point.m_z = lines[all_blocks[block][i]][j](2, 0);
                points.push_back(the_point);
            }
            insert_Line(mesh1, points, radius);
            points.clear();
        }
        mesh1.genResultMesh((".\\Results\\" + file_name + "\\" + file_block + ".obj").c_str());
        points.clear();
    }
}

void Visual::generateModelForRendering_5(vector<vector<cv::Point3d>> lines, vector<vector<cv::Point3d>> line_contain,int hegiht, int num, string file_name, int index_of_pre_node,bool judge_continue_additive, int id_continue)
{
    float radius = 0.7;
    std::vector<Vec3> points;
    string file;
    string file_2;
    string file_3;
    string file_4;
    if (judge_continue_additive == false) {
        file = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")_Layer.obj";
        file_2 = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")_Layer.gcode";
        file_3 = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")_Layer.txt";
    }
    else {
        file = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")" + "_" + to_string(id_continue) + "_Layer_subblock.obj";
        file_2 = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")" + "_" + to_string(id_continue) + "_Layer_subblock.gcode";
        file_3 = file_name + "-" + to_string(hegiht) + "_" + to_string(num) + "(" + to_string(index_of_pre_node) + ")" + "_"  + "Layer_" + to_string(id_continue)+"_subblock.txt";
    }
    ofstream file_gcode;
    file_gcode.open(file_2);
    ofstream file_txt;
    file_txt.open(file_3);
    double sum_E = 0;

    General_Mesh mesh1;
    mesh1.r = rand() / double(RAND_MAX);
    mesh1.g = rand() / double(RAND_MAX);
    mesh1.b = rand() / double(RAND_MAX);
    file_txt << lines.size() << endl;
    for (int i = 0; i < lines.size(); i++) {
        if (line_contain[i].size() == 0)
            file_txt << 0<<endl;
        else
            file_txt << 1<<endl;
        file_txt << lines[i].size() << endl;
        for (int j = 0; j < lines[i].size(); j++) {
            Vec3 the_point;
            the_point.m_x = lines[i][j].x;
            the_point.m_y = lines[i][j].y;
            the_point.m_z = lines[i][j].z;
            points.push_back(the_point);
            file_txt << the_point.m_x << " " << the_point.m_y << " " << the_point.m_z << endl;
            if (j != 0) {
                sum_E++;
                file_gcode << "G1" << " X" << the_point.m_x << " Y" << the_point.m_y << " Z" << the_point.m_z << " E" << sum_E << endl;
            }
            else
                file_gcode << "G0" << " X" << the_point.m_x << " Y" << the_point.m_y << " Z" << the_point.m_z << endl;
        }
        if (line_contain[i].size() != 0) {
            file_txt << line_contain[i].size() << endl;
            for (int j = 0; j < line_contain[i].size(); j++) {
                Vec3 the_point;
                the_point.m_x = line_contain[i][j].x;
                the_point.m_y = line_contain[i][j].y;
                the_point.m_z = line_contain[i][j].z;
                file_txt << the_point.m_x << " " << the_point.m_y << " " << the_point.m_z << endl;
            }
        }
        insert_Line(mesh1, points, radius);
        points.clear();
    }
    mesh1.genResultMesh((file).c_str());
    points.clear();

}

void Visual::generateModelForRendering_6(std::vector<cv::Point3d> orientations, string file_name)
{
    float radius = 0.1;
    General_Mesh mesh1;
    mesh1.r = 0.7;
    mesh1.g = 0.15;
    mesh1.b = 0.15;
    std::vector<Vec3> points;

    Vec3 the_point;
    Vec3 the_point_2;
    the_point.m_x = the_point.m_y = the_point.m_z = 0;
    for (int i = 0; i < orientations.size(); i++) {
        the_point.m_x += orientations[i].x;
        the_point.m_y += orientations[i].y;
        the_point.m_z += orientations[i].z;
    }
    the_point /= orientations.size();

    for (int i = 0; i < orientations.size(); i++) {
        the_point_2.m_x = orientations[i].x;
        the_point_2.m_y = orientations[i].y;
        the_point_2.m_z = orientations[i].z;

        points.push_back(the_point_2);
        points.push_back(the_point);
        insert_Line(mesh1, points, radius);
        points.clear();
    }
    mesh1.genResultMesh((file_name + "_orientation_sample_points.obj").c_str());
    points.clear();
}

void Visual::generateModelForRendering_7(cv::Point3d orientations, string file_name)
{
    float radius = 0.1;
    General_Mesh mesh1;
    mesh1.r = 0.7;
    mesh1.g = 0.15;
    mesh1.b = 0.15;
    std::vector<Vec3> points;

    Vec3 the_point;
    Vec3 the_point_2;
    the_point.m_x = the_point.m_y = the_point.m_z = 0;

    the_point_2.m_x = orientations.x * 10;
    the_point_2.m_y = orientations.y * 10;
    the_point_2.m_z = orientations.z * 10;

    points.push_back(the_point);
    points.push_back(the_point_2);
    insert_Line(mesh1, points, radius);

    mesh1.genResultMesh((file_name).c_str());
    points.clear();
}

void Visual::generateModelForRendering_8(vector<vector<Eigen::Vector3d>> lines, string file_name)
{
    float radius = 0.1;
    General_Mesh mesh1;
    mesh1.r = 0.7;
    mesh1.g = 0.15;
    mesh1.b = 0.15;
    std::vector<Vec3> points;

    Vec3 the_point;

    Vec3 the_point_2;
    for (int i = 0; i < lines.size(); i++) {
        for (int j = 0; j < lines[i].size(); j++) {
            the_point.m_x = lines[i][j].x();
            the_point.m_y = lines[i][j].y();
            the_point.m_z = lines[i][j].z();
            points.push_back(the_point);
        }
        the_point.m_x = lines[i][0].x();
        the_point.m_y = lines[i][0].y();
        the_point.m_z = lines[i][0].z();
        points.push_back(the_point);
        insert_Line(mesh1, points, radius);
        points.clear();
    }
    mesh1.genResultMesh((file_name).c_str());
}

void Visual::generateModelForRendering_9(vector<vector<cv::Point3d>> lines, string file_name)
{
    float radius = 0.3;
    General_Mesh mesh1;
    mesh1.r = 0.7;
    mesh1.g = 0.15;
    mesh1.b = 0.15;
    std::vector<Vec3> points;

    //Vec3 the_point;
    Vec3 the_point_2;
    //the_point.m_x = the_point.m_y = the_point.m_z = 0;

    for (int i = 0; i < lines.size(); i++) {
        for (int j = 0; j < lines[i].size(); j++) {
            the_point_2.m_x = lines[i][j].x;
            the_point_2.m_y = lines[i][j].y;
            the_point_2.m_z = lines[i][j].z;

            points.push_back(the_point_2);
            //points.push_back(the_point);
        }
        insert_Line(mesh1, points, radius);
        points.clear();
    }
    mesh1.genResultMesh((file_name + "_orientation_sample_points.obj").c_str());
    points.clear();
}

void Visual::generateModelForRendering_10(vector<vector<vector<Eigen::Vector3d>>> lines, vector<vector<double>>& colors, string file_name)
{
    float radius = 0.2;
    General_Mesh mesh1;
    //vector<vector<double>> colors(lines.size(), vector<double>(3));
    vector<int> num_vertex(lines.size());
    for (int t = 0; t < lines.size(); t++) {
        double r = rand() / double(RAND_MAX);
        double g = rand() / double(RAND_MAX);
        double b = rand() / double(RAND_MAX);
        
        std::vector<Vec3> points;
        Vec3 the_point;
        for (int i = 0; i < lines[t].size(); i++) {
            for (int j = 0; j < lines[t][i].size(); j++) {
                the_point.m_x = lines[t][i][j].x();
                the_point.m_y = lines[t][i][j].y();
                the_point.m_z = lines[t][i][j].z();
                points.push_back(the_point);
            }
            the_point.m_x = lines[t][i][0].x();
            the_point.m_y = lines[t][i][0].y();
            the_point.m_z = lines[t][i][0].z();
            points.push_back(the_point);
            insert_Line(mesh1, points, radius);
            
            points.clear();
        }
        num_vertex[t] = mesh1.result_vertex.size();
        colors[t][0] = r; colors[t][1] = g; colors[t][2] = b;
    }
    mesh1.genResultMesh_2(num_vertex, colors,(file_name).c_str());
}

void Visual::generateModelForRendering_11(vector<Eigen::Vector3d> points_in_cell, vector<cv::Point3d> normal, vector<vector<double>> colors, string file_name)
{
    float radius = 0.2;
    General_Mesh mesh1;
    //vector<vector<double>> colors(points_in_cell.size(), vector<double>(3));
    vector<int> num_vertex(points_in_cell.size());
    for (int t = 0; t < points_in_cell.size(); t++) {
        double r = rand() / double(RAND_MAX);
        double g = rand() / double(RAND_MAX);
        double b = rand() / double(RAND_MAX);

        std::vector<Vec3> points;
        Vec3 the_point;
        the_point.m_x = points_in_cell[t].x();
        the_point.m_y = points_in_cell[t].y();
        the_point.m_z = points_in_cell[t].z();
        points.push_back(the_point);

        Vec3 the_point_2;
        the_point_2.m_x = points_in_cell[t].x() + normal[t].x * 10;
        the_point_2.m_y = points_in_cell[t].y() + normal[t].y * 10;
        the_point_2.m_z = points_in_cell[t].z() + normal[t].z * 10;
        points.push_back(the_point_2);
        insert_Line(mesh1, points, radius);
        num_vertex[t] = mesh1.result_vertex.size();
        //colors[t][0] = r; colors[t][1] = g; colors[t][2] = b;
    }
    mesh1.genResultMesh_2(num_vertex, colors, (file_name).c_str());
}

void Visual::generateArrows(cv::Point3d normal, string file_name)
{
    float radius = 0.18;
    General_Mesh mesh1;
    //vector<vector<double>> colors(points_in_cell.size(), vector<double>(3));

    double r = rand() / double(RAND_MAX);
    double g = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);

    std::vector<Vec3> points;
    Vec3 the_point;
    the_point.m_x = 0;
    the_point.m_y = 0;
    the_point.m_z = 0;
    points.push_back(the_point);

    Vec3 the_point_2;
    the_point_2.m_x = normal.x * 1.8;
    the_point_2.m_y = normal.y * 1.8;
    the_point_2.m_z = normal.z * 1.8;
    points.push_back(the_point_2);
    insert_Line(mesh1, points, radius);

    cv::Point3d center = { the_point_2.m_x, the_point_2.m_y, the_point_2.m_z };  // 圆心
    double radius2 = 0.4;         // 半径
    int numPoints = 50;           // 生成的点数
    double translation = 0.8;

    std::vector<cv::Point3d> vertices = generatePointsOnCircle(center, normal, radius2, numPoints);

    // 平移圆心得到顶点v
    cv::Point3d  vertexV = center;
    vertexV.x += normal.x * translation;
    vertexV.y += normal.y * translation;
    vertexV.z += normal.z * translation;
    vertices.push_back(vertexV);
    vertices.push_back(center);
    // 生成圆周上点与顶点v的连接三角网格
    int vertexVIndex = numPoints;
    int centerVIndex = numPoints+1;
    std::vector<int> indices;
    for (int i = 0; i < numPoints; ++i) {
        indices.push_back(i);
        indices.push_back((i + 1) % numPoints);
        indices.push_back(vertexVIndex);

        indices.push_back(i);
        indices.push_back((i + 1) % numPoints);
        indices.push_back(centerVIndex);
    }

    saveAsObjFile(vertices, indices, file_name + "_cone.obj", mesh1, normal);

    //cout << normal.x<<"***" << endl;
   // mesh1.genResultMesh((file_name + "_array.obj").c_str());
}



void Visual::generateArrows_2(vector<Vec3>points_2, vector<cv::Point3d> normal, string file_name)
{
    float radius = 0.18;
    General_Mesh mesh1;
    //vector<vector<double>> colors(points_in_cell.size(), vector<double>(3));

    double r = rand() / double(RAND_MAX);
    double g = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);

    for (int i = 0; i < points_2.size(); i++) {
        std::vector<Vec3> points;
        Vec3 the_point;
        the_point.m_x = points_2[i].m_x;
        the_point.m_y = points_2[i].m_y;
        the_point.m_z = points_2[i].m_z;
        points.push_back(the_point);

        Vec3 the_point_2;
        the_point_2.m_x = normal[i].x * 1.8;
        the_point_2.m_y = normal[i].y * 1.8;
        the_point_2.m_z = normal[i].z * 1.8;
        points.push_back(the_point_2);
        insert_Line(mesh1, points, radius);

        cv::Point3d center = { the_point_2.m_x, the_point_2.m_y, the_point_2.m_z };  // 圆心
        double radius2 = 0.4;         // 半径
        int numPoints = 50;           // 生成的点数
        double translation = 0.8;

        std::vector<cv::Point3d> vertices = generatePointsOnCircle(center, normal[i], radius2, numPoints);

        // 平移圆心得到顶点v
        cv::Point3d  vertexV = center;
        vertexV.x += normal[i].x * translation;
        vertexV.y += normal[i].y * translation;
        vertexV.z += normal[i].z * translation;
        vertices.push_back(vertexV);
        vertices.push_back(center);
        // 生成圆周上点与顶点v的连接三角网格
        int vertexVIndex = numPoints;
        int centerVIndex = numPoints + 1;
        std::vector<int> indices;
        for (int i = 0; i < numPoints; ++i) {
            indices.push_back(i);
            indices.push_back((i + 1) % numPoints);
            indices.push_back(vertexVIndex);

            indices.push_back(i);
            indices.push_back((i + 1) % numPoints);
            indices.push_back(centerVIndex);
        }
    }
    //saveAsObjFile(vertices, indices, file_name + "_normal.obj", mesh1, normal[i]);
}