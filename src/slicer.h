#pragma once
#include<iostream>
#include<vector>
#include<array>
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include <Eigen/Dense>
using namespace std;
using VEctor = std::array<double, 3>;
using TRiangle = std::array<int, 3>;

class Slicer_2
{
public:
	// Vertex and face data.
	//using Vector = array<double, 3>;
	//using Triangle = array<int, 3>;
	std::vector<VEctor>     positions;
	std::vector<TRiangle>   triangles;
	vector<bool> jud_plane;

	// Cutting plane used for slicing.
	VEctor origin;
	VEctor normal;

	//
	// Clear mesh.
	//
	void clear()
	{
		positions.clear();
		triangles.clear();
	}

	//
	// Load OBJ file.
	//
	bool load(std::string filename)
	{
		clear();
		// Open the file.
		std::ifstream file(filename);
		if (!file) return false;
		// Read file line by line.
		std::string line;
		while (std::getline(file, line))
		{
			// Read first word in line.
			std::istringstream iss(line);
			std::string prefix;
			iss >> prefix;

			// If vertex.
			if (prefix == "v")
			{
				VEctor p;
				iss >> p[0] >> p[1] >> p[2];
				positions.push_back(p);
			}

			// If triangle.
			else if (prefix == "f")
			{
				TRiangle t;
				iss >> t[0] >> t[1] >> t[2];
				for (int i = 0; i < 3; i++)
					t[i]--;
				triangles.push_back(t);
			}

			// If something else.
			else
			{
				// We ignore the line if it's not a vertex or triangle.
			}
		}
		return true;
	}

	//
	// Save OBJ file.
	//
	bool save(std::string filename)
	{
		// Open the file.
		std::ofstream file(filename);
		if (!file) return false;
		// Write vertices.
		for (VEctor p : positions)
		{
			file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		}
		// Write triangles.
		for (TRiangle t : triangles)
		{
			for (int i = 0; i < 3; i++)
				t[i]++;
			file << "f " << t[0] << " " << t[1] << " " << t[2] << std::endl;
		}
		return true;
	}

	//
	// Slice mesh by plane.
	//
	void cut()
	{
		// We loop on all triangle indices (tid) in the mesh.
		// The do_triangle function is called on each tid,
		// and will split the triangle when required.
		//
		// If do_triangle returns false, the triangle doesn't intersect the mesh.
		// The mesh (and triangle) is unchanged and we go the next triangle (tid+=1).
		//
		// If do_triangle return true, the triangle intersects the mesh.
		// The triangle has been modified and a new triangle has been appended to the mesh.
		// Since the modified triangle might still require splitting, we stay on it (tid+=0).
		// The new appended triangle will be processed at the end.
		//
		// See do_triangle.png for a picture.

		intersections.clear();
		if (origin[0] == 0 && origin[1] == 0 && origin[2] == 0) return;
		for (size_t tid = 0; tid < triangles.size(); tid += !do_triangle(tid));
	}
	// Keep track of indices of intersection points.
	std::map<std::pair<int, int>, int> intersections;
private:

	// Floating point arithmetic is non-exact which might
	// cause problems when computing intersections.
	// We fix that using a hardcoded precision value.
	static constexpr double precision = 0.00001; //0.00001;  0.00000001

	//
	// Compute lambda such that lambda*P + (1-lambda)*Q
	// is the intersection between the plane and line PQ.
	//      lambda is in [0, 1]   =>  [PQ] intersects plane
	//      lambda is infinite    =>  [PQ] parallel to plane
	//      lambda is NaN         =>  [PQ] contained in plane
	//
	double get_lambda(VEctor p, VEctor q)
	{
		double num = (origin[0] - q[0]) * normal[0]
			+ (origin[1] - q[1]) * normal[1]
			+ (origin[2] - q[2]) * normal[2];

		double den = (p[0] - q[0]) * normal[0]
			+ (p[1] - q[1]) * normal[1]
			+ (p[2] - q[2]) * normal[2];
		if ((abs(num) < 0.000001 && abs(den) < 0.000001))//新加部分  \\if (abs(num) < 0.0000000001 && abs(den) < 0.0000000001)//新加部分
			return 0.0;
		else
			return num / den;
	}

	//
	// If edge [ij] intersects plane, add intersection vertex to mesh and return its index.
	// If the intersection vertex has already been computed before, only return its index.
	// If intersection does not exist, return -1.
	//
	int get_intersection(int i, int j)
	{
		// We require i<j.
		if (i > j) std::swap(i, j);

		// Compute lambda and return -1 if no intersection.
		VEctor p = positions[i];
		VEctor q = positions[j];
		double lambda = get_lambda(p, q);
		if (!std::isfinite(lambda) || lambda < precision || lambda > 1 - precision) return -1;

		// If the intersection has already been computed, return its index.
		auto key = std::make_pair(i, j);
		auto found = intersections.find(key);
		if (found != intersections.end()) return found->second;

		// Otherwise, compute the intersection, append it to the mesh, and return its index.
		int m = positions.size();
		positions.push_back(
			{
				lambda * p[0] + (1 - lambda) * q[0],
				lambda * p[1] + (1 - lambda) * q[1],
				lambda * p[2] + (1 - lambda) * q[2]
			});
		intersections[key] = m;
		return m;
	}

	//
	// If triangle tid intersects plane, split it and returns true.
	// Otherwise do nothing and return false.
	//
	bool do_triangle(int tid)
	{
		for (int n = 0; n < 3; n++)
		{
			int i = triangles[tid][n];
			int j = triangles[tid][(n + 1) % 3];
			int k = triangles[tid][(n + 2) % 3];

			// If edge [jk] intersects plane.
			int m = get_intersection(j, k);
			if (m != -1)
			{
				triangles.push_back({ i, j, m });
				triangles[tid] = { i, m, k };
				jud_plane.push_back(true);
				return true;
			}
		}
		return false;
	}

public:

};