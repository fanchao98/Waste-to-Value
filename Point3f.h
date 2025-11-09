#ifndef POINT3F_H
#define POINT3F_H

#include<iostream>
#include <math.h>

struct Point3ff
{
	Point3ff()
	{
		x = y = z = 0;
	}
	Point3ff(float xt, float yt, float zt)
	{
		x = xt; y = yt; z = zt;
	}
	float x, y, z;
};

struct face
{
	int v1, v2, v3;
	face(int v11 = -1, int v22 = -1, int v33 = -1)
	{
		v1 = v11;
		v2 = v22;
		v3 = v33;
	}
};


struct Point6f
{//	点的位置和向量的位置
	Point6f(float xt, float yt, float zt, float vxt, float vyt, float vzt)
	{
		position.x = xt;	position.y = yt;	position.z = zt;
		normal.x = vxt;	normal.y = vyt;	normal.z = vzt;
	}

	Point6f(Point3ff p, Point3ff n)
	{
		position.x = p.x;	position.y = p.y;	position.z = p.z;
		normal.x = n.x;		normal.y = n.y;		normal.z = n.z;
	}

	Point3ff position, normal;
};

//struct Point9f
//{//	点的位置和向量的位置
//
//	Point9f(Point3ff p, Point3ff n,Point3ff s)
//	{
//		position.x = p.x;	position.y = p.y;	position.z = p.z;
//		color.x = n.x;		color.y = n.y;		color.z = n.z;
//		size.x = s.x;	size.y = s.y;	size.z = s.z;
//	}
//
//	Point3ff position, color, size;
//};

Point3ff Cross(Point3ff A, Point3ff B);

Point3ff operator - (Point3ff A, Point3ff B);

bool operator == (Point3ff A, Point3ff B);

bool operator != (Point3ff A, Point3ff B);

Point3ff operator * (Point3ff A, float scale);

Point3ff operator + (Point3ff A, Point3ff B);

Point3ff operator / (Point3ff A, double B);

double length(Point3ff A, Point3ff B);

double length_P(Point3ff A, Point3ff B);

double Dot(Point3ff A, Point3ff B);

#endif //POINT3F_H
