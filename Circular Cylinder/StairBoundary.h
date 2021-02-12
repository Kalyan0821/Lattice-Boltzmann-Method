#ifndef STAIRBOUNDARY_H
#define STAIRBOUNDARY_H

#include <string>
#include <utility>
#include <unordered_map>
#include <set>
using namespace std;


struct Boundary{
	unordered_map<string, set<pair<int, int>>> type_to_points;
	set<pair<int, int>> boundary_points;
	set<pair<int, int>> internal_points;
};

Boundary findStairBoundary(double r, double dx, int i0, int j0);

#endif
