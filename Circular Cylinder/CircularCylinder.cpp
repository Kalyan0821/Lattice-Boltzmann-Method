#include "SolidBoundary.h"
#include<iostream>
#include<fstream>
#include<iterator>
#include<vector>
#include<cmath>
#include<numeric>
#include<algorithm>
#include<string>
using namespace std;


void writeToFile(const vector<double>& array, string name);
void writeToFile(const vector<vector<double>>& array, string name);

vector<vector<vector<double>>> collision(
	int nx, int ny, 
	const vector<vector<double>>& u, const vector<vector<double>>& v,
 	const vector<double>& cx, const vector<double>& cy, 
 	double omega, 
 	vector<vector<vector<double>>>& f, 
 	vector<vector<vector<double>>>& feq, 
 	const vector<vector<double>>& rho, 
 	const vector<double>& w,
 	int Q){

	for(int j=0; j<ny; ++j)
		for(int i=0; i<nx; ++i){
			double t1 = u[i][j]*u[i][j] + v[i][j]*v[i][j];


			for(int k=0; k<Q; ++k){
				double t2 = u[i][j]*cx[k] + v[i][j]*cy[k];

				// Lattice Boltzmann Equation- in place update				
				feq[k][i][j] = rho[i][j]*w[k]*(1. + 3.*t2 + 4.5*t2*t2 - 1.5*t1);
				f[k][i][j] = (1.-omega)*f[k][i][j] + omega*feq[k][i][j];
			}
		}
	return f;
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

vector<vector<vector<double>>> stream(vector<vector<vector<double>>>& f){
	circShift2D(f[1], 1, 0);
	circShift2D(f[2], 0, 1);
	circShift2D(f[3], -1, 0);
	circShift2D(f[4], 0, -1);
	circShift2D(f[5], 1, 1);
	circShift2D(f[6], -1, 1);
	circShift2D(f[7], -1, -1);
	circShift2D(f[8], 1, -1);
	return f;
}

vector<vector<vector<double>>> boundary(int nx, int ny, vector<vector<vector<double>>>& f, float u0, const vector<vector<double>>& rho){

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

	return f;
}


void obstacle(const Boundary& information, vector<vector<vector<double>>>& f, vector<vector<double>>& rho, string obstacle_mode){

	unordered_map<string, set<pair<int, int>>> type_to_points = information.type_to_points;
	const set<pair<int, int>>& internal_points = information.internal_points;

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

vector<vector<vector<double>>> ruv(
	int nx, int ny, 
	vector<vector<vector<double>>>& f,
	vector<vector<double>>& rho,
	vector<vector<double>>& u, vector<vector<double>>& v,
	const Boundary& information,
	int Q){

	for(int j=0; j<ny; ++j)
		for(int i=0; i<nx; ++i){

			double sum = 0.;
			for(int k=0; k<Q; ++k)
				sum += f[k][i][j];
			rho[i][j] = sum;

			u[i][j] = ((f[1][i][j]+f[5][i][j]+f[8][i][j]) - (f[3][i][j]+f[6][i][j]+f[7][i][j])) / rho[i][j]; 
			v[i][j] = ((f[2][i][j]+f[5][i][j]+f[6][i][j]) - (f[4][i][j]+f[7][i][j]+f[8][i][j])) / rho[i][j]; 
		}

	// Zero out internal points
	for(auto& p: information.internal_points){
		int i = p.first;
		int j = p.second;
		u[i][j] = 0;
		v[i][j] = 0;
		rho[i][j] = 0;
	}
	// Zero out solid points (required for bb, not required for ns)
	for(auto& p: information.solid_points){
		int i = p.first;
		int j = p.second;
		u[i][j] = 0;
		v[i][j] = 0;
	}

	return {rho, u, v};
}

int main(){

	// Choices for incompressibility, stability, boundary conditions, initialization
	double u0; // U0/c
	double alpha;  // Non-dimensional Kinematic viscosity
	int M;  // Number of cells
	int N;
	double tol;  // normalized error tolerance
	string obstacle_mode;  // set to "bb": bounce-back, "ns": no-slip (for obstacle boundary) 

	float r;  // radius, in meters
	float dx;  // cell size, in meters
	int i0;
	int j0;

	u0 = 0.15;
	alpha = 1.5;
	M = 1100;  
	N = 205;
	tol = 1e-15;
	obstacle_mode = "bb";

	r = 0.05; // meters
	dx = 1./500; // meters
	i0 = 100;
	j0 = 100;

	// Get information about solid boundary
	Boundary information = findBoundary(r, dx, i0, j0);

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

	double Dl = 2*r/dx;
	double omega = 1./(3*alpha + 0.5);  // omega = del(t)/tau

	cout << "Re: " << u0*Dl/alpha << '\n';
	cout << "Omega: " << omega << '\n';
	cout << "Ma: " << u0 * pow(3, 0.5) << '\n';
	cout << "BC: " << obstacle_mode << '\n';
	cout << "Tol: " << tol << "\n\n";

	int count = 0;
	double error = 10.;  // latest del(sum v^2)/(MN)- to determine convergence
	double ers0 = 0.;  // previous sum(v^2)/(MN)
	vector<double> errors;

	while(error>tol or count<50){

		f = collision(nx, ny, u, v, cx, cy, omega, f, feq, rho, w, Q);  // update fi values at each node using LBE
		f = stream(f);  // spatially shift fi values by ci*del(t)
		f = boundary(nx, ny, f, u0, rho);  // update values on the boundaries
		obstacle(information, f, rho, obstacle_mode);

		auto ruv_result = ruv(nx, ny, f, rho, u, v, information, Q);  // collect current results
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

		if(count==1 or count%1000==0 or (error<=tol and count >= 50))
			cout << count << ' ' << error << '\n';
	}

	// Write to file, plot with Python

	vector<double> x;  // normalized x coordinates of nodes
	vector<double> y;

	for(int i=0; i<nx; ++i)
		x.push_back(i*dx);

	for(int j=0; j<ny; ++j)
		y.push_back(j*dx);

	writeToFile(x, "./output/x.txt");
	writeToFile(y, "./output/y.txt");

	writeToFile(u, "./output/u.txt");
	writeToFile(v, "./output/v.txt");
	writeToFile(rho, "./output/rho.txt");

	writeToFile(errors, "./output/errors.txt");

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