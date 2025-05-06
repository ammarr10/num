// =======================
// File 1: LagrangeInterpolation.cpp
// =======================
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;

double lagrangeInterpolation(double x, const vector<double>& X, const vector<double>& Y) {
    int n = X.size();
    if (n != Y.size()) throw invalid_argument("X and Y must have same size.");
    
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        double term = Y[i];
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                if (X[i] == X[j]) throw invalid_argument("Duplicate X values not allowed.");
                term *= (x - X[j]) / (X[i] - X[j]);
            }
        }
        result += term;
    }
    return result;
}

/*
int main() {
    vector<double> x_vals = {1, 2, 2.5, 3, 4, 5};
    vector<double> y_vals = {1, 5, 7, 8, 2, 1};
    double interpolate_at = 3.4;
    double interp_result = lagrangeInterpolation(interpolate_at, x_vals, y_vals);
    cout << "Lagrange Interpolation at x = " << interpolate_at << " is " << interp_result << endl;
    return 0;
}
*/

// =======================
// File 2: NumericalIntegration.cpp
// =======================
#include <cmath>

double trapezoidalRule(const vector<double>& y, double h) {
    int n = y.size();
    double sum = y[0] + y[n - 1];
    for (int i = 1; i < n - 1; ++i) {
        sum += 2 * y[i];
    }
    return (h / 2.0) * sum;
}

double simpsonsRule(const vector<double>& y, double h) {
    int n = y.size();
    if ((n - 1) % 2 != 0) throw invalid_argument("Simpson's rule requires even number of intervals.");
    
    double sum = y[0] + y[n - 1];
    for (int i = 1; i < n - 1; ++i) {
        if (i % 2 == 0)
            sum += 2 * y[i];
        else
            sum += 4 * y[i];
    }
    return (h / 3.0) * sum;
}

/*
int main() {
    vector<double> fx = {-0.5, -0.4706, 0, 0.4706, 0.5, 0.2474, 0.1176};
    double a = -1, b = 2;
    int n = fx.size() - 1;
    double h = (b - a) / n;
    double trap_result = trapezoidalRule(fx, h);
    double simp_result = simpsonsRule(fx, h);
    printf("Trapezoidal Rule Result: %.6f\n", trap_result);
    printf("Simpson's 1/3 Rule Result: %.6f\n", simp_result);
    return 0;
}
*/

// =======================
// File 3: NewtonRaphson.cpp
// =======================

double f(double x) {
    return x + sin(x) - 1;
}

double f_prime(double x) {
    return 1 + cos(x);
}

double newtonRaphson(double x0, double tol_percent) {
    double x1, error;
    int max_iter = 1000;
    int iter = 0;

    do {
        double fval = f(x0);
        double deriv = f_prime(x0);
        if (fabs(deriv) < 1e-10) throw runtime_error("Derivative near zero; method may diverge.");

        x1 = x0 - fval / deriv;
        error = fabs((x1 - x0) / x1) * 100;
        x0 = x1;
        iter++;
        if (iter > max_iter) throw runtime_error("Maximum iterations exceeded.");
    } while (error > tol_percent);

    return x1;
}

/*
int main() {
    double initial_guess = 0.1;
    double tolerance = 0.5;
    double root = newtonRaphson(initial_guess, tolerance);
    cout << "Root found using Newton-Raphson: " << root << endl;
    return 0;
}
*/