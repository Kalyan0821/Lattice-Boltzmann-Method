#include "SolidBoundary.h"
#include<unordered_map>
#include<set>
#include<vector>
#include<utility>
#include<string>
#include<cmath>
#include<iostream>
#include<fstream>
using namespace std;

double rError(int i, int j, double dx, double r){

	double x = i*dx;
	double y = j*dx;
	return abs(x*x + y*y - r*r);
}

unordered_map<int, unordered_map<string, string>> getMaps(){

	unordered_map<int, unordered_map<string, string>> octant_to_type;
	octant_to_type[2] = {
		{"w", "s"}, 
		{"s", "w"}
	};
	octant_to_type[3] = {
		{"w", "s"}, 
		{"s", "e"}
	};
	octant_to_type[4] = {
		{"w", "e"}, 
		{"s", "s"}
	};
	octant_to_type[5] = {
		{"w", "e"}, 
		{"s", "n"}
	};
	octant_to_type[6] = {
		{"w", "n"}, 
		{"s", "e"}
	};
	octant_to_type[7] = {
		{"w", "n"}, 
		{"s", "w"}
	};
	octant_to_type[8] = {
		{"w", "w"}, 
		{"s", "n"}
	};

	return octant_to_type;
}

void completeBoundary(
	unordered_map<string, set<pair<int, int>>>& type_to_points,
	set<pair<int, int>>& solid_points, 
	unordered_map<int, unordered_map<string, string>>& octant_to_type){

	auto type_to_points_copy = type_to_points;

	for(auto& entry: type_to_points_copy){
		string type1 = entry.first;
		for(auto& p: type_to_points_copy[type1]){

			int i = p.first;
			int j = p.second;

			// 2
	        solid_points.insert({j, i});
	        type_to_points[octant_to_type[2][type1]].insert({j, i});

	        // 3
	        solid_points.insert({-j, i});
	        type_to_points[octant_to_type[3][type1]].insert({-j, i});

	        // 4
	        solid_points.insert({-i, j});
	        type_to_points[octant_to_type[4][type1]].insert({-i, j});

	        // 5
	        solid_points.insert({-i, -j});
	        type_to_points[octant_to_type[5][type1]].insert({-i, -j});

	        // 6
	        solid_points.insert({-j, -i});
	        type_to_points[octant_to_type[6][type1]].insert({-j, -i});

	        // 7
	        solid_points.insert({j, -i});
	        type_to_points[octant_to_type[7][type1]].insert({j, -i});

	        // 8
	        solid_points.insert({i, -j});
	        type_to_points[octant_to_type[8][type1]].insert({i, -j});
		}
	}
}


set<pair<int, int>> findInternalPoints(const set<pair<int, int>>& solid_points, int N){

	set<pair<int, int>> internal_points;
	for(int j=-(N-1); j<N; ++j)
		for(int i=0; i<N; ++i){

			if(solid_points.find({i, j}) != end(solid_points))
				break;
			else{
				internal_points.insert({i, j});
				internal_points.insert({-i, j});
			}
		}
	return internal_points;
}

set<pair<int, int>> addCenterToSet(const set<pair<int, int>>& set_pair, int i0, int j0){
	
	set<pair<int, int>> s{};

	for(auto& p: set_pair){
		int i = p.first;
		int j = p.second;
		s.insert({i+i0, j+j0});
	}

	return s;
}

void translate(
	unordered_map<string, set<pair<int, int>>>& type_to_points,
	set<pair<int, int>>& solid_points, 
	set<pair<int, int>>& internal_points,
	int i0, int j0){

	for(auto& entry: type_to_points){
		string type = entry.first;
		type_to_points[type] = addCenterToSet(type_to_points[type], i0, j0);
	}

	solid_points = addCenterToSet(solid_points, i0, j0);
	internal_points = addCenterToSet(internal_points, i0, j0);
}


Boundary findBoundary(double r, double dx, int i0, int j0){

	unordered_map<string, set<pair<int, int>>> type_to_points;
	set<pair<int, int>> solid_points;
	set<pair<int, int>> internal_points;

	int N = r/dx + 1;
	int i = N-1;
	int j = 0;
	solid_points.insert({i, j});
	type_to_points["w"].insert({i, j});
	j += 1;

	while(j <= i){
		if(rError(i, j, dx, r) <= rError(i-1, j, dx, r)){
			solid_points.insert({i, j});
        	type_to_points["w"].insert({i, j});
        }
        else if(j <= i-1 and rError(i-1, j, dx, r) < rError(i, j, dx, r)){
        	solid_points.insert({i-1, j});
	        type_to_points["w"].insert({i-1, j});
	        type_to_points["s"].insert({i, j-1});
	         
	        solid_points.insert({i-1, j-1});
	        type_to_points["w"].insert({i-1, j-1});
	        type_to_points["s"].insert({i-1, j-1});
       		i -= 1;
        }
        j += 1;
	}

	int m = round((N-1)/pow(2, 0.5));
	solid_points.insert({m, m});
	type_to_points["w"].insert({m, m});
	type_to_points["s"].insert({m, m});

	if(type_to_points["w"].find({m+1, m}) != end(type_to_points["w"]))
		type_to_points["s"].insert({m+1, m});

	unordered_map<int, unordered_map<string, string>> octant_to_type = getMaps();

	completeBoundary(type_to_points, solid_points, octant_to_type);
	internal_points = findInternalPoints(solid_points, N);

	translate(type_to_points, solid_points, internal_points, i0, j0);

	return {type_to_points, solid_points, internal_points};
}

void writeToFile(const set<pair<int, int>> solid_points, string name){
	ofstream outputfile(name);

	for(auto& p: solid_points){
		int i = p.first;
		int j = p.second;
		outputfile << i << ' ' << j << '\n';
	}
}

// int main(){

// 	Boundary information = findBoundary(0.05, 1./500, 100, 100);
// 	writeToFile(information.solid_points, "./output/solid_points.txt");
// 	writeToFile(information.internal_points, "./output/internal_points.txt");
// 	return 0;
// }

