#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#define N 30 // Maximum number of variables

typedef struct Bounds
{
    double a;
    double b;
} bounds;
struct Bounds boundlst[N];

double f(double xinputs[N]);
double integrate(double a, double n);
double montecarlo(int numvars, int iterations);
void input(int rank, int world, int *numvars, double *pointer_a, double *pointer_b, double *pointer_c, double *pointer_d, int *pointer_iterations);
void mpi_type(double *pointer_a, double *pointer_b, double *pointer_c, double *pointer_d, int *pointer_iterations, MPI_Datatype *mpi_input);
double randomgen(struct Bounds bound);
int main(int argc, char **argv)
{

    int rank, world;
    int n, local_n; // number of iterations, number of iterations per process
    int numvars;
    double local_int, total_int;
    double startwtime = 0.0, endwtime;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    input(rank, world, &numvars, &boundlst[0].a, &boundlst[0].b, &boundlst[1].a, &boundlst[1].b, &n);

    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }

    local_n = n / world; // amount of iterations each process(world) handles                                                             // locala +(number of iterations * change in x)
    local_int = montecarlo(numvars, local_n);                                             // run monte carlo
    printf("n: %d | montecarlo val: %f \n", local_n, local_int); // debug print

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

double f(double xinputs[N])
{
    // x0 = bound 1 rand num
    // x1 = bound 2 rand num
    // given 0-2, 3-7, 4-5 from inside out evaluate to ~197.177
    return pow(xinputs[0], 3) + sin(xinputs[1]) + pow(2, xinputs[3]);
}

double montecarlo(int numvars, int iterations)
{

    double fVal;
    double sum = 0;
    int fVals[iterations];

    int curIt = 0;

    while (curIt < iterations - 1)
    {
        // Sample the function's values at the random point
        double xinputs[30];
        for(int i = 0; i<numvars; i++){
            xinputs[i] = randomgen(boundlst[i]);
        }
        fVal = f(xinputs);
        fVals[curIt] = fVal;
        // Add the f(x) value to the running sum
        sum += fVal;

        curIt++;
    }
    double mean = (sum / iterations);
    float stnderr = 0;
    for (int i = 0; i < iterations; i++)
    {
        stnderr += pow(fVals[i] - mean, 2);
    }
    stnderr = sqrt((stnderr / iterations)) / sqrt(iterations);

    printf("Standard Error: %f\n", stnderr);
    double estimate = mean;
    // (b-a)*(d-c)*...*mean = integral
    for(int i = 0; i < numvars; i++){
        estimate *= (boundlst[i].b-boundlst[i].a);
    }
    return estimate;
}

void mpi_type(double *pointer_a, double *pointer_b, double *pointer_a1, double *pointer_b1, int *pointer_iterations, MPI_Datatype *mpi_input)
{
    MPI_Aint address_a, address_b, address_a1, address_b1, addresss_iterations;
    MPI_Aint displacement_arr[5] = {0};
    int blocklengths_arr[5] = {1, 1, 1, 1, 1};
    MPI_Datatype types_arr[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
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

void input(int rank, int world, int *numvars, double *pointer_a, double *pointer_b, double *pointer_a1, double *pointer_b1, int *pointer_iterations)
{
    MPI_Datatype mpi_input;
    mpi_type(pointer_a, pointer_b, pointer_b1, pointer_b1, pointer_iterations, &mpi_input);

    if (rank == 0)
    {
        printf("Enter number of variables: \n");
        scanf("%d", numvars);
        printf("Enter number of iterations: \n");
        scanf("%d", pointer_iterations);
        for (int i = 0; i < *numvars; i++)
        {
            printf("Enter bounds for bound %d: ", i);
            scanf("%lf %lf", &boundlst[i].a, &boundlst[i].b);
        }
    }

    MPI_Bcast(pointer_a, 1, mpi_input, 0, MPI_COMM_WORLD);

    MPI_Type_free(&mpi_input);
}

