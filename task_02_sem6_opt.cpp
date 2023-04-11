#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <mpi.h>

using namespace std;

void gaussian_elimination(double *arr, double *part_arr, int n, int rank, int size) {
    int i, j, k;
    double factor;
    MPI_Status status;

    for (k = 0; k < n - 1; k++) {
        int root = (k * size) / n;

        if (rank == root) {
            for (i = 0; i < n; i++) {
                arr[k * n + i] = part_arr[(k % (n / size)) * n + i];
            }
        }

        MPI_Bcast(arr + k * n, n, MPI_DOUBLE, root, MPI_COMM_WORLD);

        for (i = k + 1; i < n; i++) {
            if ((i * size) / n == rank) {
                factor = part_arr[(i % (n / size)) * n + k] / arr[k * n + k];
                for (j = k + 1; j < n; j++) {
                    part_arr[(i % (n / size)) * n + j] -= factor * arr[k * n + j];
                }
                part_arr[(i % (n / size)) * n + k] = 0;
            }
        }
    }
}
int get_index_p(int g_line,int size){
    return g_line/(size-1);
}
int get_thread_num(int g_line,int size){
    int p_line=(g_line+1)%(size-1);
    if (p_line==0){
        p_line=size-1;
    }
    return p_line;
}
void printmtx(double *arr, int32_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << arr[i * n + j] << " ";
        }
        cout << endl;
    }
}

void printpartmtx(double *arr, int row, int n){
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << arr[i*n + j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Wrong input" << endl;
        return -1;
    }
    string input_file = argv[1];
    ifstream mtx(input_file);
    if (!mtx.is_open()) {
        cerr << "Error opening file: " << input_file << endl;
        return 1;
    }

    int n;
    mtx >> n;
    int g_line=0;
    int rank = 0, size = 0;
    double start_time = 0, end_time = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;
    double *arr = new double[n * n];
    double *part_arr = new double[n*(n/(size-1))];

    if (rank == 0) {
        for (int i = 0; i < n * n; ++i) {
            mtx >> arr[i];
        }

        cout << "threads: "<< size << endl;
        cout << n << endl;
        printmtx(arr, n);
        cout << endl;
//send
    int thread_index=1;
        while (g_line<n){
            if (thread_index == size){
                thread_index=1;
            }
            MPI_Send(&arr[g_line*n], n, MPI_DOUBLE, thread_index, get_index_p(g_line,size), comm);
            thread_index++;
            g_line++;
        }
    }
    if (rank){
        for (int i = 0; i < n/(size-1); ++i) {
            MPI_Recv(&part_arr[i*n],n,MPI_DOUBLE,0,i,comm,&status);
        }
    }
    MPI_Barrier(comm);
    if(rank == 3){
        printpartmtx(part_arr, size-1, n);
        cout << endl;
    }
    g_line=0;
    if (rank==0){

    } else{


    }


//    MPI_Scatter(arr, n, MPI_DOUBLE, part_arr, n, MPI_DOUBLE, 0, comm);
//    start_time = MPI_Wtime();
//    gaussian_elimination(arr, part_arr, n, rank, size);
//    end_time = MPI_Wtime();
//    MPI_Gather(part_arr, n, MPI_DOUBLE, arr, n, MPI_DOUBLE, 0, comm);
//
//    if (rank == 0) {
//        cout << "Upper triangular matrix:" << endl;
//        printmtx(arr, n);
//        cout << "Calculation time: " << end_time - start_time << "s" << endl;
//    }

    delete[] arr;
    delete[] part_arr;
    mtx.close();

    MPI_Finalize();
    return 0;
}






