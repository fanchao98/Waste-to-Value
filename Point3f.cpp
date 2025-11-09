#include"Point3f.h"


Point3ff Cross(Point3ff A, Point3ff B)
{
	return Point3ff(A.y * B.z - A.z * B.y, B.x * A.z - A.x * B.z, A.x * B.y - A.y * B.x);
}

//A-B
Point3ff operator - (Point3ff A, Point3ff B)
{
	return Point3ff(A.x - B.x, A.y - B.y, A.z - B.z);
}


bool operator == (Point3ff A, Point3ff B)
{
	return (A.x == B.x) && (A.y == B.y) && (A.z == B.z);
}

bool operator !=(Point3ff A, Point3ff B)
{
	return (A.x != B.x) || (A.y != B.y) || (A.z != B.z);
}

Point3ff operator * (Point3ff A, float scale)
{
	return Point3ff(A.x * scale, A.y * scale, A.z * scale);
}

Point3ff operator + (Point3ff A, Point3ff B)
{
	return Point3ff(A.x + B.x, A.y + B.y, A.z + B.z);
}

double length(const Point3ff A, const Point3ff B)
{
	return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z);
}


double length_P(const Point3ff A, const Point3ff B)
{
	return sqrtf((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z));
}

Point3ff operator / (Point3ff A, double B) {
	return Point3ff(A.x / B, A.y / B, A.z / B);
}

double Dot(Point3ff A, Point3ff B)
{
	return A.x * B.x + A.y * B.y + A.z * B.z;
}


