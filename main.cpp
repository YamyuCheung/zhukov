#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <mpi.h>

using namespace std;

void gaussian_elimination(double* part_arr, int n, int rank, int size) {
    double factor;
    int global, str_rank, str_num;
    for (global = 0; global < n; ++global) {
        str_rank = global % (size - 1) + 1;
        str_num = global / (size - 1);

        if (rank == str_rank) {
            MPI_Bcast(&part_arr[str_num * n], n, MPI_DOUBLE, str_rank, MPI_COMM_WORLD);
        }

        for (int i = 0; i < n / (size - 1); ++i) {
            if (i != str_num) {
                factor = part_arr[i * n + global] / part_arr[str_num * n + global];
                for (int j = global + 1; j < n; ++j) {
                    part_arr[i * n + global] -= factor * part_arr[str_num * n + global];
                }
            }
        }
    }

}
void printmtx(double* arr, int32_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << arr[i * n + j] << " ";
        }
        cout << endl;
    }
}

void printpartmtx(double* arr, int row, int n) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << arr[i * n + j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {
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
    int rank = 0, size = 0, dest = 1;
    double start_time = 0, end_time = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;
    double* arr = new double[n * n];
    double* part_arr = new double[n * (n / (size - 1))];
    double* line_arr = new double[n];

    if (rank == 0) {
        for (int i = 0; i < n * n; ++i) {
            mtx >> arr[i];
        }
        for (int i = 0; i < n; ++i) {
            if (dest == size) {
                dest = 1;
            }
            MPI_Send(arr + n * i, n, MPI_DOUBLE, dest, 0, comm);
            dest++;
        }

        cout << "processes: " << size << endl;
        cout << n << endl;
        printmtx(arr, n);
        cout << endl;
    }
    if (rank != 0) {
        for (int j = 0; j < n / (size - 1); ++j) {
            MPI_Recv(part_arr + n * j, n, MPI_DOUBLE, 0, 0, comm, &status);
        }
    }
    MPI_Barrier(comm);
    if (rank == 3) {
        printpartmtx(part_arr, size - 1, n);
        cout << endl;
    }
    if (rank != 0) {
        gaussian_elimination(part_arr, n, rank, size);
        MPI_Barrier(comm);
        //        for (int i = 0; i < n; ++i) {
        //            if(rank == i%(size-1)+1){
        //                MPI_Send(&part_arr[(i/(size-1))*n+0], n, MPI_DOUBLE, 0, i, comm);
        //            }
        //        }
        //        MPI_Barrier(comm);
        if (rank == 3) {
            printpartmtx(part_arr, size - 1, n);
            cout << endl;
        }
    }
    else {
        //        for (int i = 0; i < n; ++i) {
        //            MPI_Recv(&arr[i*n+0], n, MPI_DOUBLE, i%3+1, i, comm, &status);
        //        }

    }
    //    printmtx(arr, n);


    delete[] arr;
    delete[] part_arr;
    delete[] line_arr;
    mtx.close();

    MPI_Finalize();
    return 0;
}






