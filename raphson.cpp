#include <iostream>
#include <cmath>
using namespace std;

double f(double x) {
    return x + sin(x) - 1;
}

double f_prime(double x) {
    return 1 + cos(x);
}

double newtonRaphson(double x0, double tol_percent) {
    double x1, error;
    do {
        x1 = x0 - f(x0) / f_prime(x0);
        error = fabs((x1 - x0) / x1) * 100;
        x0 = x1;
    } while (error > tol_percent);
    return x1;
}

int main() {
    double initial_guess = 0.1;
    double tolerance = 0.5;
    double root = newtonRaphson(initial_guess, tolerance);
    cout << "The root of the equation is " << root;
    return 0;
}