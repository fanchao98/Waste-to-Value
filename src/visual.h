#pragma once
#ifndef VISUAL_H_
#define VISUAL_H_

#include "config.h"
#include <vector>

#include <exception>
#include <fstream>

#include <Eigen/Dense>
#include "slicer.h"
#include"GeneralMesh.h"
#include<string>
#include <time.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;
class Visual
{
public:
    struct Circle {
        cv::Point3d center;
        cv::Point3d normal;
        double radius;
    };

    void insert_Line(General_Mesh& mesh, std::vector<Vec3> points, float radius);
    void generateModelForRendering(std::vector<std::vector<std::vector<Eigen::MatrixXd>>> lines, string file_name);
    void generateModelForRendering_2(Eigen::Vector3d lines, string file_name);
    void generateModelForRendering_3(Eigen::Vector3d lines, string file_name, std::vector<Eigen::MatrixXd> vis_points);
    void generateModelForRendering_4(std::vector<std::vector<Eigen::MatrixXd>> lines, string file_name, vector<vector<int>> all_blocks);
    void generateModelForRendering_5(vector<vector<cv::Point3d>> lines, vector<vector<cv::Point3d>> line_contain, int hegiht,int num, string file_name, int index_of_pre_node, bool judge_continue_additive, int id_continue);
    void generateModelForRendering_6(std::vector<cv::Point3d> orientations, string file_name);
    void generateModelForRendering_7(cv::Point3d orientations, string file_name);
    void generateModelForRendering_8(vector<vector<Eigen::Vector3d>> lines, string file_name);
    void generateModelForRendering_9(vector<vector<cv::Point3d>> lines, string file_name);
    void generateModelForRendering_10(vector<vector<vector<Eigen::Vector3d>>> lines, vector<vector<double>>& colors, string file_name);
    void generateModelForRendering_11(vector<Eigen::Vector3d> points, vector<cv::Point3d> normal, vector<vector<double>> colors,string file_name);
    void generateArrows(cv::Point3d normal, string file_name);
    void generateArrows_2(vector<Vec3>points_2, vector<cv::Point3d> normal, string file_name);

    std::vector<cv::Point3d> generatePointsOnCircle(const cv::Point3d center, const cv::Point3d normal, double radius, int numPoints) {
        std::vector<cv::Point3d> points;
        double angleIncrement = 2 * 3.14159 / numPoints;

        for (int i = 0; i < numPoints; i++) {
            double angle = i * angleIncrement;
            double x = radius * cos(angle);
            double y = radius * sin(angle);
            points.push_back({ x, y, 0 });
        }

        std::vector<Eigen::MatrixXd> temp_V;
        temp_V.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
            temp_V[i].resize(3, 1);
            temp_V[i](0, 0) = points[i].x;
            temp_V[i](1, 0) = points[i].y;
            temp_V[i](2, 0) = points[i].z;
        }
        Eigen::Vector3d vectorAfter(0, 0, 1);
        Eigen::Vector3d vectorBefore(normal.x, normal.y, normal.z);
        Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
        for (int i = 0; i < numPoints; i++)
            temp_V[i] = rotMatrix.inverse() * temp_V[i];
        for (int i = 0; i < numPoints; i++) {
            points[i].x = temp_V[i](0, 0) + center.x;
            points[i].y = temp_V[i](1, 0) + center.y;
            points[i].z = temp_V[i](2, 0) + center.z;
        }

        return points;
    };

    void saveAsObjFile(const std::vector<cv::Point3d>& vertices, const std::vector<int>& indices, const std::string filename, General_Mesh mesh1, cv::Point3d normal) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Failed to create file: " << filename << std::endl;
            return;
        }

        for (const cv::Point3d& vertex : vertices) {
            file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
        }

        for (int i = 0; i < indices.size(); i += 3) {
            file << "f " << indices[i]+1 << " " << indices[i + 1] + 1 << " " << indices[i + 2] + 1 << std::endl;
        }

        for (int i = 0; i < mesh1.result_vertex.size(); i++)
            file << "v " << mesh1.result_vertex[i].m_x << " " << mesh1.result_vertex[i].m_y << " " << mesh1.result_vertex[i].m_z << endl;
        for (int i = 0; i < mesh1.result_faces.size(); i++)
        {
            int facet_count = mesh1.result_faces[i].size();
            file << "f";
            for (int j = 0; j < facet_count; j++)
            {
                file << " " << mesh1.result_faces[i][j] + 1 + vertices.size();
            }
            file << endl;
        }
        file <<"%%% " << normal.x << " " << normal.y << " " << normal.z << endl;
    }


    
};
#endif