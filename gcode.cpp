#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <array>
#include <math.h>
#include <exception>
#include <fstream>

#include "datastructures.h"
#include "config.h"
#include "katana.h"
#include "stl.h"
#include "gcode.h"
#include "polygon.h"
using namespace std;
std::vector<std::vector<bool>> is_flatten_area;
std::vector<pair<int, int>> index_flatten_layer;

bool check_inside(Point pt, Point* pgn_begin, Point* pgn_end, K traits)
{
    switch (CGAL::bounded_side_2(pgn_begin, pgn_end, pt, traits)) {
    case CGAL::ON_BOUNDED_SIDE:
        return true;
        break;
    case CGAL::ON_BOUNDARY:
        return false;
        break;
    case CGAL::ON_UNBOUNDED_SIDE:
        return false;
        break;
    }
}

void adjust_gcode(vector<vector<vector<Vertex>>>& ALL_OPP)
{
    Vertex center_point;
    center_point.x = center_point.y = center_point.z = 0;
    int cont_point = 0;
    for (int i = 0; i < ALL_OPP.size(); i++)
        for (int j = 0; j < ALL_OPP[i].size(); j++)
            for (int k = 0; k < ALL_OPP[i][j].size(); k++)
            {
                center_point.x += ALL_OPP[i][j][k].x;
                center_point.y += ALL_OPP[i][j][k].y;
                center_point.z += ALL_OPP[i][j][k].z;
                cont_point++;
            }
    center_point.x /= cont_point;
    center_point.y /= cont_point;
    center_point.z /= cont_point;
    for (int i = 0; i < ALL_OPP.size(); i++)
        for (int j = 0; j < ALL_OPP[i].size(); j++)
            for (int k = 0; k < ALL_OPP[i][j].size(); k++)
            {
                ALL_OPP[i][j][k].x -= center_point.x;
                ALL_OPP[i][j][k].y -= center_point.y;
                ALL_OPP[i][j][k].z -= center_point.z;
            }
    for (int i = 0; i < ALL_OPP.size(); i++)
        for (int j = 0; j < ALL_OPP[i].size(); j++)
            for (int k = 0; k < ALL_OPP[i][j].size(); k++)
            {
                double temp = ALL_OPP[i][j][k].y;
                ALL_OPP[i][j][k].y = ALL_OPP[i][j][k].z;
                ALL_OPP[i][j][k].z = -temp;
            }
    for (int i = 0; i < ALL_OPP.size(); i++)
        for (int j = 0; j < ALL_OPP[i].size(); j++)
            for (int k = 0; k < ALL_OPP[i][j].size(); k++)
            {
                ALL_OPP[i][j][k].x += 100;
                ALL_OPP[i][j][k].y += 100;
                ALL_OPP[i][j][k].z += 6.5;
            }
}

double distance(Vertex a, Vertex b)
{
    double dis = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
    return dis;
}

// save Gcode
// iterates over the previously generated layers and emit gcode for every segment
// uses some configuration values to decide when to retract the filament, how much
// to extrude and so on.
//void GCode::write(const char* filename, std::vector<Layer>& layers, float min_z)
void GCodeWriter::write(vector<vector<Vertex>>& all_slice_points)
{
    all_slice_points.clear();

    // segments shorter than this are ignored
    float skipDistance = .0001;

    // offset of the emitted Gcode coordinates to the .stl ones
    //Vertex offset={75,75,Katana::Instance().config.get("z_offset")-Katana::Instance().min_z};
    Vertex offset = { 0,0,0 };
    Vertex position = { 0,0,0 };
    vector<vector<Vertex>> normal_segments;
    
    normal_segments.resize(Katana::Instance().layers.size());

    for (unsigned int i = 0; i < Katana::Instance().layers.size(); i++)
    {
        vector<Vertex> temp;
        //all_slice_points.push_back(temp);
        Layer& l = Katana::Instance().layers[i];

        float feedrate = (i == 0) ? 500.f : 1200.f;
        //fprintf(file, "G1 Z%f F%f ;layer %d\n", l.z + offset.z, feedrate, i); // move to layer's z plane

        float extrusion = (i == 0) ? 1 : 0; // extrusion axis position

        normal_segments[i].resize(l.segments.size());
        for (unsigned int j = 0; j < l.segments.size(); j++)
        {
            normal_segments[i][j] = l.segments[j].normal;
            Vertex& v0 = l.segments[j].vertices[0];
            Vertex& v1 = l.segments[j].vertices[1];

            //assert(v0.z==l.z);
            //assert(v1.z==l.z);
            // skip segment with NaN or Infinity caused by numeric instablities
            if (v0.z != l.z) continue;
            if (v1.z != l.z) continue;

            if (v1.distance(position) < v0.distance(position))
                std::swap(v0, v1);


            if (j != l.segments.size() - 1)
            {
                if (v0 == l.segments[j + 1].vertices[0] || v0 == l.segments[j + 1].vertices[1])
                    swap(v0, v1);
            }


            // check distance to decide if we need to travel
            float d = v0.distance(position);
            if (d > skipDistance)
            {
                vector<Vertex> temp2;
                bool temp_bool = true;
                all_slice_points.push_back(temp2);

                Vertex VV;
                VV.x = v0.x + offset.x; VV.y = v0.y + offset.y; VV.z = l.z + offset.z;
                all_slice_points[all_slice_points.size() - 1].push_back(VV);

                position = v0;
            }

            if (v1.distance(position) > skipDistance)
            {
                Vertex VV;
                VV.x = v1.x + offset.x; VV.y = v1.y + offset.y; VV.z = l.z + offset.z;
                all_slice_points[all_slice_points.size() - 1].push_back(VV);

                position = v1;
            }
        }

        if (all_slice_points[all_slice_points.size() - 1].size() == 0) {
            all_slice_points.erase(all_slice_points.end() - 1);
        }
    }
    ////////////////////////////////////////////////////////////////////////
}
