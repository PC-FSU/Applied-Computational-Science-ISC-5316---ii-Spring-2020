#include<mpi.h>
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>

#define XMAX 20.0
#define XMIN -20.0
#define YMAX 20.0
#define YMIN -20.0
#define CLOUDSIZE 1000
#define DIMENSION 2
#define ITERATIONS 10000

double random_double(double min, double max);
int random_int(int min, int max);

using namespace std;

int main(int argc, char* argv[]) {

    // MPI initialization
    int npes,mype,ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    // test function
    auto f = [](double x, double y) -> double {
        return 1.0 + pow(sin(x), 2) + pow(sin(y), 2) - 0.1 * exp(-pow(x, 2) -pow(y, 2));
    };

    // initialize cloud of randomly chosen points on a single thread
    double Vx[CLOUDSIZE];
    double Vy[CLOUDSIZE];
    double Vf[CLOUDSIZE];
    if(mype == 0) {
        for(int i = 0; i < CLOUDSIZE; i++) {
            Vx[i] = random_double(XMIN, XMAX);
            Vy[i] = random_double(YMIN, YMAX);
            Vf[i] = f(Vx[i], Vy[i]);
        }
    }

    for(int l = 0; l < ITERATIONS; l+=npes) {
        
        // bcast cloud from root thread to all others -- repeat at each iteration since multiple points may have changed
        MPI_Bcast(Vx, CLOUDSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(Vy, CLOUDSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(Vf, CLOUDSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // find maximum value f(M)
        int max_index = 0;
        for(int i = 1; i < CLOUDSIZE; i++) {
            if (Vf[i] > Vf[max_index]) {
                max_index = i;
            }
        }

        // choose n+1 points from the cloud at random
        int simplex[DIMENSION + 1];
        for(int i = 0; i < DIMENSION + 1; i++) {
            simplex[i] = random_int(0, CLOUDSIZE - 1);
        }

        // sort on function value
        sort(
            simplex,
            simplex + DIMENSION + 1,
            [&Vf](const int a, const int b) {
                return Vf[a] > Vf[b];
            }
        );

        // find the centroid G
        double Gx;
        double Gy;
        double Gf;
        double sum = 0;
        for(int i = 0; i < DIMENSION; i++) {
            sum += Vx[simplex[i]];
        }
        Gx = sum / DIMENSION;
        sum = 0;
        for(int i = 0; i < DIMENSION; i++) {
            sum += Vy[simplex[i]];
        }
        Gy = sum / DIMENSION;
        Gf = f(Gx, Gy);

        // project from the maximum point of the simplex through the centroid
        double Px;
        double Py;
        double Pf;
        Px = 2 * Gx - Vx[simplex[DIMENSION]];
        Py = 2 * Gy - Vy[simplex[DIMENSION]];
        Pf = f(Px, Py);

        // keep if the point is less than the global maximum
        if(XMIN <= Px && Px <= XMAX && YMIN <= Py && Py <= YMAX) {
            if(Pf < Vf[max_index]) {
                Vx[max_index] = Px;
                Vy[max_index] = Py;
                Vf[max_index] = Pf;
                // cout << "Thread " << mype << " sending Pf = " << Pf << endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // send updates
        if(mype == 0) {
            for(int i = 1; i < npes; i++) {
                double temp_Vx, temp_Vy, temp_Vf;
                MPI_Recv(&temp_Vx, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&temp_Vy, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&temp_Vf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // double-check success condition prior to updating root's cloud
                for(int i = 1; i < CLOUDSIZE; i++) {
                    if (Vf[i] > Vf[max_index]) {
                        max_index = i;
                    }
                }
                if(temp_Vf < Vf[max_index]) {
                    // cout << "Replacing max f " << Vf[max_index] << " by " << temp_Vf << " at index " << max_index << endl;
                    Vx[max_index] = temp_Vx;
                    Vy[max_index] = temp_Vy;
                    Vf[max_index] = temp_Vf;
                }
            }
        }
        else {
            MPI_Send(&Vx[max_index], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&Vy[max_index], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&Vf[max_index], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

    }

    MPI_Finalize();

    // print resulting cloud w/ function values
    if(mype == 0) {
        for(int i = 0; i < CLOUDSIZE; i++) {
            cout << Vf[i] << "\t" << Vx[i] << "\t" << Vy[i] << endl;
        }
    }

}

// generate random double in range
double random_double(double min, double max) {

    // dude seriously fuck this
    mt19937 generator(chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> dist(min, max);
    return dist(generator);

}

// generate random integer in range
int random_int(int min, int max) {

    // get your shit together, C++
    mt19937 generator(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> dist(min, max);
    return dist(generator);

}
