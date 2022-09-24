#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#define N 2 // Number of variables

typedef struct Bounds
{
    double a;
    double b;
} bounds;
struct Bounds boundlst[N];

double f(double x0, double x1);
double integrate(double a, double n);
double montecarlo(double a, double b, int iterations);
void input(int rank, int world, double *pointer_a, double *pointer_b, double *pointer_c, double *pointer_d, int *pointer_iterations);
void mpi_type(double *pointer_a, double *pointer_b, double *pointer_c, double *pointer_d, int *pointer_iterations, MPI_Datatype *mpi_input);
double randomgen(struct Bounds bound);
int main(int argc, char **argv)
{

    int rank, world;
    int n, local_n; // number of iterations, number of iterations per process
    double local_int, total_int;
    double a, b, c, d, dx, local_a, local_b; // start and end of integral             // change in x
    double startwtime = 0.0, endwtime;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    input(rank, world, &boundlst[0].a, &boundlst[0].b, &boundlst[1].a, &boundlst[1].b, &n);

    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }

    dx = (b - a) / n;
    local_n = n / world; // amount of iterations each process(world) handles

    local_a = a + rank * local_n * dx;                                                             // start of integral + (rank# * number of iterations * change in x)
    local_b = local_a + local_n * dx;                                                              // locala +(number of iterations * change in x)
    local_int = montecarlo(local_a, local_b, local_n);                                             // run monte carlo
    printf("n: %d | a: %f | b: %f | montecarlo val: %f \n", local_n, local_a, local_b, local_int); // debug print

    MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // end timer
        endwtime = MPI_Wtime();
        printf("wall clock: %lf\n", endwtime - startwtime);
    }
    MPI_Finalize();
}

double randomgen(struct Bounds bound)
{

    float ran = (float)(rand());
    double val = bound.a + (ran / RAND_MAX) * (bound.b - bound.a);

    return val;
}

double f(double x0, double x1)
{
    // x0 = bound 1 rand num
    // x1 = bound 2 rand num
    return pow(x0, 3) + sin(x1);
}

double montecarlo(double a, double b, int iterations)
{

    double fVal;
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
    }

    int curIt = 0;

    while (curIt < iterations - 1)
    {
        // Sample the function's values at the random point
        fVal = f(randomgen(boundlst[0]), randomgen(boundlst[1]));

        // Add the f(x) value to the running sum
        sum += fVal;

        curIt++;
    }

    double estimate = (b - a) * (sum / iterations);

    return estimate;
}

void mpi_type(double *pointer_a, double *pointer_b, double *pointer_a1, double *pointer_b1, int *pointer_iterations, MPI_Datatype *mpi_input)
{
    MPI_Aint address_a, address_b, address_a1, address_b1, addresss_iterations;
    MPI_Aint displacement_arr[5] = {0};
    int blocklengths_arr[5] = {1, 1, 1, 1, 1};
    MPI_Datatype types_arr[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Get_address(pointer_a, &address_a);
    MPI_Get_address(pointer_b, &address_b);
    MPI_Get_address(pointer_a1, &address_a1);
    MPI_Get_address(pointer_b1, &address_b1);
    MPI_Get_address(pointer_iterations, &addresss_iterations);
    displacement_arr[1] = address_b - address_a;
    displacement_arr[2] = address_a1 - address_a;
    displacement_arr[3] = address_b1 - address_a;
    displacement_arr[4] = addresss_iterations - address_a;
    MPI_Type_create_struct(5, blocklengths_arr, displacement_arr, types_arr, mpi_input);
    MPI_Type_commit(mpi_input);

    // displacement_arr[2] = addresss_iterations - address_a;
    // MPI_Type_create_struct(3, blocklengths_arr, displacement_arr, types_arr, mpi_input);
    // MPI_Type_commit(mpi_input);
}

void input(int rank, int world, double *pointer_a, double *pointer_b, double *pointer_a1, double *pointer_b1, int *pointer_iterations)
{
    MPI_Datatype mpi_input;
    mpi_type(pointer_a, pointer_b, pointer_b1, pointer_b1, pointer_iterations, &mpi_input);

    if (rank == 0)
    {
        printf("Enter number of iterations\n");
        scanf("%d", pointer_iterations);
        for (int i = 0; i < N; i++)
        {
            printf("Enter bounds for bound %d: ", i);
            scanf("%lf %lf", &boundlst[i].a, &boundlst[i].b);
        }
    }

    MPI_Bcast(pointer_a, 1, mpi_input, 0, MPI_COMM_WORLD);

    MPI_Type_free(&mpi_input);
}