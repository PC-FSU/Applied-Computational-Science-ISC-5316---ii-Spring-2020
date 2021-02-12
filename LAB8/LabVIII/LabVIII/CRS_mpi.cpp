#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "mpi.h"

using std::cout;
using std::endl;
using std::pow;
using std::sin;
using std::exp;
using std::max_element;
using std::min_element;
using std::ofstream;
using std::calloc;


double f(double x, double y) {
	return 1 + pow(sin(x), 2) + pow(sin(y), 2) - 0.1*exp(-x*x - y*y);
}

double* generateCloud(int N, double lower, double upper) {
	unsigned seed = 10;//std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(lower, upper);
	
	double* cloud = new double[N*2];
	for (int i = 0; i < 2*N; i++) {
		cloud[i] = distribution(generator);
	}
	return cloud;
}

int main(int argc, char *argv[]) {
	//initialize MPI
	int id, Np;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &Np);
	
	//main program variables
	int N = 1000;
	int n = 2;
	int max_iter = 10000;
	double* cloud = (double*)calloc(N*2, sizeof(double));
	double* fcloud = (double*)calloc(N, sizeof(double));
	double fM_prev, fM_curr;
	int Midx;
	int idx;
	double max_vertex_f, max_vertex_x, max_vertex_y;
	double accumx, accumy;
	double centroidx, centroidy;
	double P[3];
	double global_min;
	int i_suc = 0;
	bool converg = false;
	double diff;
	bool new_max = false;
	bool new_P = false;
	
	//variables for parallelization
	int pts_per_proc;
	double* buffer = (double*)calloc(N, sizeof(double));
	//double* P_ext = (double*)calloc(3, sizeof(double));
	double P_ext[3]={0};
	MPI_Status status;
	
	//random integer generator
	unsigned seed = 10;//std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<int> distribution(0, N-1);
	

	ofstream fout("max-min.dat");

	auto start = std::chrono::high_resolution_clock::now();
     
	//(master) create cloud and broadcast
	if (id == 0) {
		cloud = generateCloud(N, -20, 20);
	}
	MPI_Bcast(cloud, N*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//(all) evaluate objective function at assigned points and combine results
	pts_per_proc = N / Np;
	for (int i = 0; i < pts_per_proc; i++) {
		buffer[id*pts_per_proc + i] = f(cloud[n*id*pts_per_proc + i*n], cloud[n*id*pts_per_proc + i*n + 1]);
	}
	MPI_Allreduce(buffer, fcloud, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	
	fM_curr = *max_element(fcloud, fcloud + N);

	//main loop
	while (!converg && i_suc <= max_iter) {
		//find centroid and "maximum" vertex
		accumx = 0;
		accumy = 0;
		for (int j = 0; j < n+1; j++) {
			idx = distribution(generator);
			if (j == 0) {
				max_vertex_f = fcloud[idx];
				max_vertex_x = cloud[idx*n];
				max_vertex_y = cloud[idx*n + 1];
			}
			else {
				if (fcloud[idx] > max_vertex_f) {
					accumx += max_vertex_x;
					accumy += max_vertex_y;
					max_vertex_f = fcloud[idx];
					max_vertex_x = cloud[idx*n];
					max_vertex_y = cloud[idx*n + 1];
				}
				else {
					accumx += cloud[idx*n];
					accumy += cloud[idx*n + 1];
				}
			}
		}
		
		centroidx = accumx / n;
		centroidy = accumy / n;
		//find P
		P[0] = 2*centroidx - max_vertex_x;
		P[1]= 2*centroidy - max_vertex_y;
		P[2] = f(P[0], P[1]);
		
		//find M
		for (int j = 0; j < N; j++) {
			if (fcloud[j] == fM_curr) {
				Midx = j;
				break;
			}
		}
		
		//check internal P
		//if valid -> swap points and send to other process
		if ((P[0] > -20 && P[0] < 20) && (P[1] > -20 && P[1] < 20) && (P[2] < fM_curr)) {
			cloud[Midx*n] = P[0];
			cloud[Midx*n + 1] = P[1];
			fcloud[Midx] = P[2];

			new_P = true;
			for (int j = 0; j < Np; j++) {
				if (j != id) {
					MPI_Send(&new_P, 1, MPI_CXX_BOOL, j, 100, MPI_COMM_WORLD);
					MPI_Send(P, 3, MPI_DOUBLE, j, 200, MPI_COMM_WORLD);
				}
			}
			/*
			//difference between max and min
			if (i_suc == 0) {
				fM_curr = *max_element(fcloud, fcloud + N);
			}
			global_min = *min_element(fcloud, fcloud + N);
			fout << fM_curr - global_min << endl;
			*/
			new_max = true;
			i_suc++;
		}
		else {
			new_P = false;
			for (int j = 0; j < Np; j++) {
				if (j != id) {
					cout<<"here";
					MPI_Send(&new_P, 1, MPI_CXX_BOOL, j, 100, MPI_COMM_WORLD);
				}
			}
		}
		
		//receive external P -> check duplicate -> swap
		for (int j = 0; j < Np; j++) {
			if (j != id) {
				MPI_Recv(&new_P, 1, MPI_CXX_BOOL, j, 100, MPI_COMM_WORLD, &status);
				
				if (new_P) {
					//receive external P
					MPI_Recv(P_ext, 3, MPI_DOUBLE, j, 200, MPI_COMM_WORLD, &status);
					//find M
					fM_curr = *max_element(fcloud, fcloud + N);
					for (int k = 0; k < N; k++) {
						if (fcloud[k] == fM_curr) {
							Midx = k;
							break;
						}
					}
					//swap
					cloud[Midx*n] = P_ext[0];
					cloud[Midx*n + 1] = P_ext[1];
					fcloud[Midx] = P_ext[2];
					/*
					//difference between max and min
					global_min = *min_element(fcloud, fcloud + N);
					fout << fM_curr - global_min << endl;
					*/
					new_max = true;
					i_suc++;
				}
			}
		}
		
		if (new_max) {
			fM_prev = fM_curr;
			fM_curr = *max_element(fcloud, fcloud + N);
			diff = fabs(fM_curr - fM_prev)/fM_curr;
			
			if (diff < pow(10, -16) && diff > 0) {
				converg = true;
				MPI_Bcast(&converg, 1, MPI_CXX_BOOL, id, MPI_COMM_WORLD);
			}
			
			new_max = false;
		}
		
		//cout << id << ':' << Midx << endl;
	}

	auto stop = std::chrono::high_resolution_clock::now();

	auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	
	global_min = *min_element(fcloud, fcloud + N);
	cout << runtime.count() << '|' << global_min << endl;

	
	fout.close();
	MPI_Finalize();
	
}

