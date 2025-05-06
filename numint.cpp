#include <iostream>
#include <vector>
using namespace std;

double trapezoidalRule(const vector<double> y, double h) {
    int n = y.size();
    double sum = y[0] + y[n - 1];
    for (int i = 1; i < n - 1; ++i) {
        sum += 2 * y[i];
    }
    return (h / 2.0) * sum;
}

double simpsonsRule(const vector<double> y, double h) {
    int n = y.size();
    double sum = y[0] + y[n - 1];
    for (int i = 1; i < n - 1; i++) {
        if (i % 2 == 0)
            sum += 2 * y[i];
        else
            sum += 4 * y[i];
    }
    return (h / 3.0) * sum;
}

int main() {
    vector<double> fx = {-0.5, -0.4706, 0, 0.4706, 0.5, 0.2474, 0.1176};
    double a = -1, b = 2;
    int n = fx.size() - 1;
    double h = (b - a) / n;
    double trap = trapezoidalRule(fx, h);
    double simp = simpsonsRule(fx, h);
    printf("Trapezoidal Rule Result: %.6f\n", trap);
    printf("Simpson's 1/3 Rule Result: %.6f\n", simp);
    return 0;
}