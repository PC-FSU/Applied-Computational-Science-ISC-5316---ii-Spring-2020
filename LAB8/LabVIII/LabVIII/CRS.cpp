#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>

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
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(lower, upper);
	double* cloud = new double[N*2];
	for (int i = 0; i < 2*N; i++) {
		cloud[i] = distribution(generator);
	}
	return cloud;
}

int main() {
	int N = 1000;
	int n = 2;
	int max_iter = 10000;
	double* cloud;
	double fcloud[N];
	double fM_prev, fM_curr;
	int Midx;
	int idx;
	double max_vertex_f, max_vertex_x, max_vertex_y;
	double accumx, accumy;
	double centroidx, centroidy;
	double fP, Px, Py;
	double global_min;
	int i_suc = 0;
	bool converg = false;
	double diff;
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<int> distribution(0, N-1);
	
	ofstream fout("cloud.dat");
	
	auto start = std::chrono::high_resolution_clock::now();

	//create cloud and evaluate objective function at all points
	cloud = generateCloud(N, -20, 20);
	for (int i = 0; i < N; i++) {
		fcloud[i] = f(cloud[i*n], cloud[i*n + 1]);
	}
	
	
	fM_curr = *max_element(fcloud, fcloud + N);
	
	//main loop
	while (!converg && i_suc <= max_iter) {
		//store cloud at particular iterations
		if (i_suc == 999 || i_suc == 4999 || i_suc == 7499 || i_suc == 9999) {
			for (int j = 0; j < N; j++) {
				fout << cloud[j*n] << ' ';
			}
			fout << endl;
			for (int j = 0; j < N; j++) {
				fout << cloud[j*n + 1] << ' ';
			}
			fout << endl;
		}
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
		Px = 2*centroidx - max_vertex_x;
		Py = 2*centroidy - max_vertex_y;
		fP = f(Px, Py);
		
		//find M
		for (int j = 0; j < N; j++) {
			if (fcloud[j] == fM_curr) {
				Midx = j;
				break;
			}
		}
		
		//check P and maybe update
		if ((Px > -20 && Px < 20) && (Py > -20 && Py < 20) && (fP < fM_curr)) {
			fcloud[Midx] = fP;
			cloud[Midx*n] = Px;
			cloud[Midx*n + 1] = Py;

			fM_prev = fM_curr;
			fM_curr = *max_element(fcloud, fcloud + N);
			diff = fabs(fM_curr - fM_prev)/fM_curr;
			if (diff < pow(10, -18) && diff > 0) {
				converg = true;
			}
			i_suc++;
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();

	auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	
	global_min = *min_element(fcloud, fcloud + N);
	cout << runtime.count() << '|' << global_min << endl;
	
	fout.close();
}
