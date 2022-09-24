#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 2
typedef struct Bounds
{
    double a;
    double b;
} bounds;
struct Bounds boundlst[N];

double f(double x0, double x1);
double montecarlo(double a, double b, int iterations);
double randomgen(struct Bounds bound);
void input();

int main(int argc, char argv)
{
    input();
    double result = montecarlo(boundlst[0].a, boundlst[0].b, 1000000);
    printf("Result: %f \t", result);

    return 0;
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
    // time_t t;
    // srand((unsigned)time(&t));
    double random, fVal;
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
    }

    int curIt = 0;

    while (curIt < iterations - 1)
    {

        // Select a random number within the limits of integration
        float ran = (float)(rand());
        random = a + (ran / RAND_MAX) * (b - a);

        // Sample the function's values
        // double poiinta1 = boundlst[0].a;
        double b1 = boundlst[0].b;
        int a2 = boundlst[1].a;
        int b2 = boundlst[1].b;
        fVal = f(randomgen(boundlst[0]), randomgen(boundlst[1]));

        // Add the f(x) value to the running sum
        sum += fVal;

        curIt++;
    }

    double estimate = (b - a) * (sum / iterations);

    return estimate;
}

void input()
{
    int iterations;
    printf("Enter number of iterations:");
    scanf("%d", &iterations);
    for (int i = 0; i < N; i++)
    {
        printf("Enter bounds for bound %d: ", i);
        scanf("%lf %lf", &boundlst[i].a, &boundlst[i].b);
    }

    // printf("List of bounds: ");
    // printf("bound 0 %lf %lf", boundlst[0].a, boundlst[0].b);
    // printf("bound 1 %lf %lf", boundlst[1].a, boundlst[1].b);
    // printf("Iterations: %d", iterations);
}