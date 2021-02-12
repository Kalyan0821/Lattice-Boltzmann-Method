#ifndef CURVEDBOUNDARY_H
#define CURVEDBOUNDARY_H

#include <vector>
#include <utility>

using namespace std;


class SurroundingPoint{
public:
	pair<int, int> pointF;  // Point in fluid/on boundary (i, j)
	pair<int, int> pointB;  // Connected point inside solid/on boundary
	pair<int, int> pointFF; // Point in fluid on F-B line beyond F
	double delta;           // Fraction of F-B distance inside fluid, will be in [0, 1]
	int dir_i;              // Direction into solid
	int dir_j;				// Opposite to i

	SurroundingPoint(int iF, int jF, int iB, int jB, int N, int dir_i);
	void setDelta(int N);
	void addCenter(int i0, int j0);  // Call to translate circle to required origin
};

vector<SurroundingPoint> getSurroundingPoints(double r, double dx, int i0, int j0, bool add_center=true);
double solveQuadratic(double a, double b, double c);
int oppositeDirection(int i);

#endif
