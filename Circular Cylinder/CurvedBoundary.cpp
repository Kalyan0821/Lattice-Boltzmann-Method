#include "CurvedBoundary.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


SurroundingPoint::SurroundingPoint(int iF, int jF, int iB, int jB, int N, int dir_i){
	pointF = {iF, jF};
	pointB = {iB, jB};
	pointFF = {2*iF-iB, 2*jF-jB};  // F is the midpoint of B and FF
	setDelta(N);
	this->dir_i = dir_i;
	this->dir_j = oppositeDirection(dir_i);
}

void SurroundingPoint::setDelta(int N){
	int iF = pointF.first;
	int jF = pointF.second;
	int iB = pointB.first;
	int jB = pointB.second;

	// Set up quadratic equation to find delta
	double a = (iB-iF)*(iB-iF) + (jB-jF)*(jB-jF);
	double b = 2 * (iF*(iB-iF) + jF*(jB-jF));
	double c = iF*iF + jF*jF - (N-1)*(N-1);

	this->delta = solveQuadratic(a, b, c);
}

void SurroundingPoint::addCenter(int i0, int j0){
	pointF.first += i0;
	pointF.second += j0;

	pointB.first += i0;
	pointB.second += j0;

	pointFF.first += i0;
	pointFF.second += j0;
}


double solveQuadratic(double a, double b, double c){
	double D = b*b - 4*a*c;
	double delta1 = (-b + sqrt(D))/(2*a);
	double delta2 = (-b - sqrt(D))/(2*a);

	// Return the solution in [0, 1]
	if(0<=delta1 and delta1<=1)
		return delta1;
	else if(0<=delta2 and delta2<=1)
		return delta2;
}

int oppositeDirection(int i){
	if(i==1)
		return 3;
	if(i==2)
		return 4;
	if(i==3)
		return 1;
	if(i==4)
		return 2;
	if(i==5)
		return 7;
	if(i==6)
		return 8;
	if(i==7)
		return 5;
	if(i==8)
		return 6;
}

bool isSolid(int i, int j, int N){
	return (i*i + j*j <= (N-1)*(N-1));
}

bool isFluid(int i, int j, int N){
	return (i*i + j*j >= (N-1)*(N-1));
}

vector<vector<int>> getSolidNeighbors(int i, int j, int N){

	vector<vector<int>> candidates = {{4, i, j-1},  // dir_i, (coordinates)
									  {2, i, j+1},
									  {3, i-1, j},
									  {7, i-1, j-1},
									  {6, i-1, j+1},
									  {1, i+1, j},
									  {8, i+1, j-1},
									  {5, i+1, j+1}};
	
	vector<vector<int>> solid_neighbors;
	for(auto& p: candidates){
		int dir_i = p[0];
		int i = p[1];
		int j = p[2];

		if(isSolid(i, j, N))
			solid_neighbors.push_back({dir_i, i, j});
	}

	return solid_neighbors;
}


vector<SurroundingPoint> getSurroundingPoints(double r, double dx, int i0, int j0, bool add_center/*=true*/){
	int N = r/dx + 1;  // Number of points in [(0, 0), (R, 0)]

	vector<SurroundingPoint> surrounding_points;
	for(int i=-N; i<=N; ++i)
		for(int j=-N; j<=N; ++j){

			if(isFluid(i, j, N)){

				vector<vector<int>> solid_neighbors = getSolidNeighbors(i, j, N);
				
				for(auto& p: solid_neighbors){
					int dir_i = p[0];
					int iB = p[1];
					int jB = p[2];
					
					SurroundingPoint surrounding_point(i, j, iB, jB, N, dir_i);
					if(add_center)
						surrounding_point.addCenter(i0, j0);
					surrounding_points.push_back(surrounding_point);
				}
			}
		}

	return surrounding_points;
}

void writePointsToFile(const vector<SurroundingPoint>& surrounding_points, string nameF, string nameB, string nameFF, string name_delta){
	ofstream outF(nameF);
	ofstream outB(nameB);
	ofstream outFF(nameFF);
	ofstream out_delta(name_delta);


	for(auto& p: surrounding_points){
		outF << p.pointF.first << ' ' << p.pointF.second << '\n';
		outB << p.pointB.first << ' ' << p.pointB.second << '\n';
		outFF << p.pointFF.first << ' ' << p.pointFF.second << '\n';
		out_delta << p.delta << '\n';
	}
}


// int main(){
// 	vector<SurroundingPoint> surrounding_points = getSurroundingPoints(0.05, 1./1000, 200, 200);
// 	writePointsToFile(surrounding_points, "./notebooks/F.txt", "./notebooks/B.txt", "./notebooks/FF.txt", "./notebooks/delta.txt");
// 	return 0;
// }

