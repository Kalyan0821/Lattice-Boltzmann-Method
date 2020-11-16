#ifndef SOLIDBOUNDARY_H
#define SOLIDBOUNDARY_H

#include<string>
#include<unordered_map>
#include<set>
using namespace std;


struct Boundary{
	unordered_map<string, set<pair<int, int>>> type_to_points;
	set<pair<int, int>> solid_points;
	set<pair<int, int>> internal_points;
};

Boundary findBoundary(double r, double dx, int i0, int j0);


double radiusError(int i, int j, double dx, double r);

unordered_map<int, unordered_map<string, string>> getMaps();

void completeBoundary(
	unordered_map<string, set<pair<int, int>>>& type_to_points,
	set<pair<int, int>>& solid_points, 
	unordered_map<int, unordered_map<string, string>>& octant_to_type);

set<pair<int, int>> findInternalPoints(const set<pair<int, int>>& solid_points, int N);

set<pair<int, int>> addCenterToSet(const set<pair<int, int>>& set_pair, int i0, int j0);

void translate(
	unordered_map<string, set<pair<int, int>>>& type_to_points,
	set<pair<int, int>>& solid_points, 
	set<pair<int, int>>& internal_points,
	int i0, int j0);

#endif
