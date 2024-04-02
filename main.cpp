#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

vector<double> f(vector<double> x) {
    return {8 * x[0] * x[0] + 4 * x[0] * x[1] + 5 * x[1] * x[1] - 2 * x[1],
            x[0] - x[1]};
}

vector<double> gradient(vector<double> x) {
    return {16 * x[0] + 4 * x[1], 4 * x[0] + 10 * x[1] - 2};
}

vector<vector<double>> hessian() {
    return {{16, 4}, {4, 10}};
}

vector<double> solve_linear_system(vector<vector<double>> hess, vector<double> grad) {
    int n = hess.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double factor = hess[j][i] / hess[i][i];
            for (int k = i; k < n; k++) {
                hess[j][k] -= factor * hess[i][k];
            }
            grad[j] -= factor * grad[i];
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += hess[i][j] * x[j];
        }
        x[i] = (grad[i] - sum) / hess[i][i];
    }

    return x;
}

// Newton's method implementation
vector<double> newton(vector<double> x, double eps) {
    vector<double> grad = gradient(x);
    vector<vector<double>> hess = hessian();

    while (sqrt(grad[0] * grad[0] + grad[1] * grad[1]) >= eps) {

        cout << "x" << 0 << " = [" << x[0] << ", " << x[1] << "]^T" << endl <<
             "delf" " = [" << grad[0] << ", " << grad[1] << "]" << endl <<
             "||delf||" " = " << sqrt(grad[0]*grad[0] + grad[1]*grad[1]) << endl;

        vector<double> d = solve_linear_system(hess, grad);

        for (int i = 0; i < x.size(); i++) {
            x[i] -= d[i];
        }

        cout <<  "x" << 1 << " = [" << x[0] << ", " << x[1] << "]" << endl <<
                 "f" << 1 << " = " << f(x)[0] << endl;

        grad = gradient(x);
    }

    return x;
}

int main() {
    vector<double> x0 = {6, 6};
    double eps = 1e-6;

    vector<double> solution = newton(x0, eps);

    cout << "\nMinimum point: x = (" << fixed << setprecision(3) << solution[0] << ", " << solution[1] << "), Minimum value: " << f(solution)[0] << endl;

    cout << endl;

    return 0;
}
