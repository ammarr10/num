#include <iostream>
#include <vector>
using namespace std;

double lagrangeInterpolation(double x, vector<double> X, vector<double> Y) {
    int n = X.size();
    for (int i = 0; i < n; ++i) {
        if (X[i] > x) {
            double x0 = X[i - 1], x1 = X[i];
            double y0 = Y[i - 1], y1 = Y[i];
            double term1 = y0 * ((x - x1) / (x0 - x1));
            double term2 = y1 * ((x - x0) / (x1 - x0));
            return term1 + term2;
        }
    }
    throw invalid_argument("x is out of range of X values.");
}

int main() {
    vector<double> x = {1, 2, 2.5, 3, 4, 5};
    vector<double> y = {1, 5, 7, 8, 2, 1};
    double value = 3.4;
    double interpolated = lagrangeInterpolation(value, x, y);
    cout << "The interpolated value at x = " << value << " is " << interpolated;
    return 0;
}