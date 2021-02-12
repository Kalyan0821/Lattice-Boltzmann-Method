// To compile and run:
// g++ -std=c++17 -fopenmp CircularCylinder.cpp StairBoundary.cpp CurvedBoundary.cpp -o outP
// ./outP


#include "StairBoundary.h"
#include "CurvedBoundary.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <string>
#include <omp.h>
using namespace std;

#define PARALLELIZE

void writeToFile(const vector<double>& array, string name);
void writeToFile(const vector<vector<double>>& array, string name);

void collision(
	int nx, int ny, 
	const vector<vector<double>>& u, const vector<vector<double>>& v,
 	const vector<double>& cx, const vector<double>& cy, 
 	double omega, 
 	vector<vector<vector<double>>>& f, 
 	vector<vector<vector<double>>>& feq, 
 	const vector<vector<double>>& rho, 
 	const vector<double>& w,
 	int Q){

	#ifdef PARALLELIZE
	#pragma omp parallel for collapse(2)
	#endif
	for(int j=0; j<ny; ++j)
		for(int i=0; i<nx; ++i){
			double t2 = u[i][j]*u[i][j] + v[i][j]*v[i][j];


			for(int k=0; k<Q; ++k){
				double t1 = cx[k]*u[i][j] + cy[k]*v[i][j];

				// Lattice Boltzmann Equation- in place update				
				feq[k][i][j] = rho[i][j]*w[k]*(1. + 3.*t1 - 1.5*t2 + 4.5*t1*t1);
				f[k][i][j] = (1.-omega)*f[k][i][j] + omega*feq[k][i][j];
			}
		}
}


void circShift2D(vector<vector<double>>& a, int s0, int s1){

	if(s0>0)
		rotate(begin(a), end(a)-s0, end(a));  // Shift rows by s0 steps down
	else
		rotate(begin(a), begin(a)-s0, end(a));

	for(auto& r: a){
		if(s1>0)
			rotate(begin(r), end(r)-s1, end(r));  // Shift columns by s1 steps to the right
		else
			rotate(begin(r), begin(r)-s1, end(r));
	}
}

void stream(vector<vector<vector<double>>>& f){
	circShift2D(f[1], 1, 0);
	circShift2D(f[2], 0, 1);
	circShift2D(f[3], -1, 0);
	circShift2D(f[4], 0, -1);
	circShift2D(f[5], 1, 1);
	circShift2D(f[6], -1, 1);
	circShift2D(f[7], -1, -1);
	circShift2D(f[8], 1, -1);
}

void boundary(int nx, int ny, vector<vector<vector<double>>>& f, float u0, const vector<vector<double>>& rho){

	int N = ny-1;

	// West boundary
	for(int j=1; j<ny-1; ++j){

		double uj = 6*u0*(j)*(N-j)/(N*N);

		f[1][0][j] = f[3][0][j] + 2.*rho[0][j] * uj/3.;

		f[5][0][j] = f[7][0][j] - 0.5*(f[2][0][j]-f[4][0][j]) + rho[0][j] * uj/6.;
		f[8][0][j] = f[6][0][j] + 0.5*(f[2][0][j]-f[4][0][j]) + rho[0][j] * uj/6.;
	}
	
	// East boundary
	for(int j=0; j<ny; ++j){

		f[3][nx-1][j] = f[3][nx-2][j];
		f[7][nx-1][j] = f[7][nx-2][j];
		f[6][nx-1][j] = f[6][nx-2][j];
	}

	// South boundary
	for(int i=0; i<nx; ++i){

		f[2][i][0] = f[4][i][0];
		f[5][i][0] = f[7][i][0];
		f[6][i][0] = f[8][i][0];
	}

	// North boundary
	for(int i=0; i<nx; ++i){

		f[4][i][ny-1] = f[2][i][ny-1];
		f[7][i][ny-1] = f[5][i][ny-1];
		f[8][i][ny-1] = f[6][i][ny-1];
	}
}


void stairBoundaryUpdate(const Boundary& information, vector<vector<vector<double>>>& f, string obstacle_mode){

	unordered_map<string, set<pair<int, int>>> type_to_points = information.type_to_points;

	for(auto& p: type_to_points["w"]){
		int i = p.first;
		int j = p.second;

		if(obstacle_mode=="bb"){
			f[1][i][j] = f[3][i][j];
			f[5][i][j] = f[7][i][j];
			f[8][i][j] = f[6][i][j];
		}
		else if(obstacle_mode=="ns"){
			f[1][i][j] = f[3][i][j];
			f[5][i][j] = f[7][i][j] - 0.5*(f[2][i][j] - f[4][i][j]);
			f[8][i][j] = f[6][i][j] + 0.5*(f[2][i][j] - f[4][i][j]);
		}
	}

	for(auto& p: type_to_points["e"]){
		int i = p.first;
		int j = p.second;

		if(obstacle_mode=="bb"){
			f[3][i][j] = f[1][i][j];
			f[7][i][j] = f[5][i][j];
			f[6][i][j] = f[8][i][j];
		}
		else if(obstacle_mode=="ns"){
			f[3][i][j] = f[1][i][j];
			f[7][i][j] = f[5][i][j] + 0.5*(f[2][i][j]-f[4][i][j]);
			f[6][i][j] = f[8][i][j] - 0.5*(f[2][i][j]-f[4][i][j]);			
		}
	}

	for(auto& p: type_to_points["s"]){
		int i = p.first;
		int j = p.second;

		if(obstacle_mode=="bb"){
			f[2][i][j] = f[4][i][j];
			f[5][i][j] = f[7][i][j];
			f[6][i][j] = f[8][i][j];
		}
		else if(obstacle_mode=="ns"){
			f[2][i][j] = f[4][i][j];
			f[5][i][j] = f[7][i][j] - 0.5*(f[1][i][j]-f[3][i][j]);
			f[6][i][j] = f[8][i][j] + 0.5*(f[1][i][j]-f[3][i][j]);			
		}
	}

	for(auto& p: type_to_points["n"]){
		int i = p.first;
		int j = p.second;

		if(obstacle_mode=="bb"){
			f[4][i][j] = f[2][i][j];
			f[7][i][j] = f[5][i][j];
			f[8][i][j] = f[6][i][j];
		}
		else if(obstacle_mode=="ns"){
			f[4][i][j] = f[2][i][j];
			f[7][i][j] = f[5][i][j] + 0.5*(f[1][i][j]-f[3][i][j]);
			f[8][i][j] = f[6][i][j] - 0.5*(f[1][i][j]-f[3][i][j]);		
		}
	}

}

void curvedBoundaryUpdate(
	const vector<SurroundingPoint>& surrounding_points,
	vector<vector<vector<double>>>& f,
    const vector<vector<double>>& u, const vector<vector<double>>& v,
 	const vector<double>& cx, const vector<double>& cy,
 	const vector<vector<double>>& rho,
 	const vector<double>& w,
 	double omega,
 	int curved_scheme){

	for(auto& surrounding_point: surrounding_points){

		int iF = surrounding_point.pointF.first;
		int jF = surrounding_point.pointF.second;
		int iB = surrounding_point.pointB.first;
		int jB = surrounding_point.pointB.second;
		int iFF = surrounding_point.pointFF.first;
		int jFF = surrounding_point.pointFF.second;

		double delta = surrounding_point.delta;
		int dir_i = surrounding_point.dir_i;
		int dir_j = surrounding_point.dir_j;

		double X;
		double uBF, vBF;

		if(delta >= 0.5){
			uBF = (1-1./delta)*u[iF][jF] + (1./delta)*0;  // Zero wall velocity
			vBF = (1-1./delta)*v[iF][jF] + (1./delta)*0;
			X = omega * (2*delta-1);
		}
		else{

			if(curved_scheme==1){
				uBF = u[iF][jF];
				vBF = v[iF][jF];
				X = (omega/(1-omega)) * (2*delta-1);
			}

			else if(curved_scheme==2){
				uBF = u[iFF][jFF];
				vBF = v[iFF][jFF];
				X = (omega/(1-2*omega)) * (2*delta-1);
			}
		}

		double t1 = cx[dir_i]*uBF + cy[dir_i]*vBF;
		double t2 = u[iF][jF]*u[iF][jF] + v[iF][jF]*v[iF][jF];
		double t3 = cx[dir_i]*u[iF][jF] + cy[dir_i]*v[iF][jF];

		double fi_star =  rho[iF][jF]*w[dir_i]*(1. + 3.*t1 - 1.5*t2 + 4.5*t3*t3);  // Zero wall velocity

		f[dir_j][iB][jB] = (1-X)*f[dir_i][iF][jF] + (X)*fi_star - 0; 
	}
}

vector<vector<vector<double>>> ruv(
	int nx, int ny, 
	vector<vector<vector<double>>>& f,
	vector<vector<double>>& rho,
	vector<vector<double>>& u, vector<vector<double>>& v,
	const Boundary& information,
	string obstacle_mode,
	int Q){

	#ifdef PARALLELIZE
	#pragma omp parallel for collapse(2)
	#endif
	for(int j=0; j<ny; ++j)
		for(int i=0; i<nx; ++i){

			double sum = 0.;
			for(int k=0; k<Q; ++k)
				sum += f[k][i][j];
			rho[i][j] = sum;

			u[i][j] = ((f[1][i][j]+f[5][i][j]+f[8][i][j]) - (f[3][i][j]+f[6][i][j]+f[7][i][j])) / rho[i][j]; 
			v[i][j] = ((f[2][i][j]+f[5][i][j]+f[6][i][j]) - (f[4][i][j]+f[7][i][j]+f[8][i][j])) / rho[i][j]; 
		}

	if(obstacle_mode=="curved"){
		// Zero out internal points
		for(auto& p: information.internal_points){
			int i = p.first;
			int j = p.second;
			u[i][j] = 0;
			v[i][j] = 0;
			rho[i][j] = 0;
		}
	}

	else{
		// Zero out internal points
		for(auto& p: information.internal_points){
			int i = p.first;
			int j = p.second;
			u[i][j] = 0;
			v[i][j] = 0;
			rho[i][j] = 0;
		}
		// Zero out boundary points (not required for ns, but required for bb since it does not zero out tangential-to-surface components)
		for(auto& p: information.boundary_points){
			int i = p.first;
			int j = p.second;
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}

	return {rho, u, v};
}

int main(){

	// Given specifications (in meters)
	double H = 0.41;  
	double L = 2.2;   
	double D = 0.1;   
	double x0 = 0.2;  
	double y0 = 0.2;  

	// Incompressibility, stability, convergence
	double u0;   // U0/c
	int N0;      // D/dx, diameter in lattice units
	double nul;  // non-dimensional kinematic viscosity
	double tol;  // normalized error tolerance

	// Obstacle boundary conditions
	string obstacle_mode;  // set to "bb": bounce-back, "ns": no-slip (for obstacle boundary), "curved": curved
	int curved_scheme;  // set to 1 or 2 for "curved" obstacle modes

	// Choices
	u0 = 0.15;
	N0 = 50;
	nul = 1.5;
	obstacle_mode = "curved";
	curved_scheme = 1;
	
	tol = 1e-15;  


	double r = D/2.;  // radius, in meters
	double dx = D/N0;  // cell size, in meters
	int M = L/dx;  // Number of cells
	int N = H/dx;
	int i0 = x0/dx;
	int j0 = y0/dx;
	double omega = 1./(3*nul + 0.5);  // omega = del(t)/tau


	Boundary information = findStairBoundary(r, dx, i0, j0);
	vector<SurroundingPoint> surrounding_points = getSurroundingPoints(r, dx, i0, j0);

	//////////////////////////////////////////////////////////////////////////////

	// D2Q9 arrangement
	int Q = 9; 
	vector<double> w{4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};  // weighting factors
	vector<double> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};  // velocities
	vector<double> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
	double cs2 = 1./3;

	int nx = M+1;  // # nodes = (1 + # lattices) 
	int ny = N+1;  


	vector<vector<double>> u(nx, vector<double>(ny, 0));
	vector<vector<double>> v(nx, vector<double>(ny, 0));
	vector<vector<double>> rho(nx, vector<double>(ny, 1)); 
	vector<vector<vector<double>>> f(Q, vector<vector<double>>(nx, vector<double>(ny)));
	vector<vector<vector<double>>> feq(Q, vector<vector<double>>(nx, vector<double>(ny)));

	cout << "Re: " << u0*N0/nul << '\n';
	cout << "Nu: " << nul <<", Omega: " << omega << '\n';
	cout << "Ma: " << u0 * pow(3, 0.5) << '\n';
	cout << "BC: " << obstacle_mode << '\n';
	if(obstacle_mode=="curved")
		cout << "Curved Scheme: " << curved_scheme << '\n';
	cout << "Tol: " << tol << "\n\n";

	int count = 0;
	double error = 10.;  // latest del(sum v^2)/(MN)- to determine convergence
	double ers0 = 0.;  // previous sum(v^2)/(MN)
	vector<double> errors;

	while(error>tol or count<50){

		if(obstacle_mode=="curved"){
			collision(nx, ny, u, v, cx, cy, omega, f, feq, rho, w, Q);
			curvedBoundaryUpdate(surrounding_points, f, u, v, cx, cy, rho, w, omega, curved_scheme);
			stream(f);
			boundary(nx, ny, f, u0, rho);
		}

		else{
			collision(nx, ny, u, v, cx, cy, omega, f, feq, rho, w, Q);  // update fi values at each node using LBE
			stream(f);  // spatially shift fi values by ci*del(t)
			boundary(nx, ny, f, u0, rho);  // update values on the boundaries
			stairBoundaryUpdate(information, f, obstacle_mode);
		}


		auto ruv_result = ruv(nx, ny, f, rho, u, v, information, obstacle_mode, Q);  // collect current results
		rho = ruv_result[0];
		u = ruv_result[1];
		v = ruv_result[2];

		++count;

		double ers = 0.;  // current sum(v^2)/(MN)
		for(int i=0; i<nx; ++i)
			for(int j=0; j<ny; ++j)
				ers += u[i][j]*u[i][j] + v[i][j]*v[i][j];

		ers /= (M*N);

		error = abs(ers-ers0);
		ers0 = ers;

		errors.push_back(error);

		if(count==1 or count%1==0 or (error<=tol and count >= 50))
			cout << count << ' ' << error << '\n';
	}

	// Write to file, plot with Python

	vector<double> x;  // normalized x coordinates of nodes
	vector<double> y;

	for(int i=0; i<nx; ++i)
		x.push_back(i*dx);

	for(int j=0; j<ny; ++j)
		y.push_back(j*dx);


	string s = "";
	if(obstacle_mode=="curved")
		s = "_curved";
	writeToFile(x, "./output"+s+"/x.txt");
	writeToFile(y, "./output"+s+"/y.txt");
	writeToFile(u, "./output"+s+"/u.txt");
	writeToFile(v, "./output"+s+"/v.txt");
	writeToFile(rho, "./output"+s+"/rho.txt");
	writeToFile(errors, "./output"+s+"/errors.txt");

	return 0;
}

void writeToFile(const vector<double>& array, string name){
	ofstream outputfile(name);

	for(auto e: array)
		outputfile << e << ' ';
}

void writeToFile(const vector<vector<double>>& array, string name){
	ofstream outputfile(name);

	for(auto r: array){
		for(auto e: r)
			outputfile << e << ' ';

		outputfile << '\n';
	}
}


