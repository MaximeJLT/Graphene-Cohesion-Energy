#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

const double PI = 3.14159265358979323846;

double cosd(double angle) {
    return std::cos(angle * PI / 180.0);
}

double sind(double angle) {
    return std::sin(angle * PI / 180.0);
}

int main() {
    double pm = 2.47;
    int N = 15;

    std::vector<double> a = {pm, 0.0, 0.0};
    std::vector<double> b = {pm * cosd(60), pm * sind(60), 0.0};
    std::vector<double> c = {0.0, 0.0, 1.0};

    std::string name = "C";
    int Na = 2 * (N + 1) * (N + 1);
    std::vector<std::vector<double>> R_n(Na, std::vector<double>(3, 0.0));

    double d = pm / (2.0 * cosd(30)); // interatomic distance
    double D = std::sqrt(3.0 / 4.0) * d; // hole-bridge distance
    double A = 15.2; // eV
    double B = 24100;

    std::vector<double> motif_1 = {0.0, 0.0, 0.0};
    std::vector<double> motif_2 = {pm / 2.0, pm / 2.0 * std::tan(PI / 6.0), 0.0};

    int n = 0;
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (n < Na) {
                R_n[n][0] = i * a[0] + j * b[0] + motif_1[0];
                R_n[n][1] = i * a[1] + j * b[1] + motif_1[1];
                R_n[n][2] = i * a[2] + j * b[2] + motif_1[2];
                n++;
            }
            if (n < Na) {
                R_n[n][0] = i * a[0] + j * b[0] + motif_2[0];
                R_n[n][1] = i * a[1] + j * b[1] + motif_2[1];
                R_n[n][2] = i * a[2] + j * b[2] + motif_2[2];
                n++;
            }
        }
    }

    // Extract X_s and Y_s
    std::vector<double> X_s, Y_s;
    for (const auto& atom : R_n) {
        X_s.push_back(atom[0]);
        Y_s.push_back(atom[1]);
    }

    // Determine adsorption site coordinates
    double X_H = *max_element(X_s.begin(), X_s.end()) / 2.0;
    double Y_H = *max_element(Y_s.begin(), Y_s.end()) / 2.0;

    double X_T = X_H;
    double Y_T = Y_H + d;

    double X_B = X_H + D;
    double Y_B = Y_H;

    // Generate XYZ file
    std::ofstream file("Graphene.xyz");
    if (file.is_open()) {
        file << Na << "\nGraphene\n";
        for (const auto& atom : R_n) {
            file << std::setw(3) << name << " "
                 << std::fixed << std::setprecision(3)
                 << std::setw(8) << atom[0] << " "
                 << std::setw(8) << atom[1] << " "
                 << std::setw(8) << atom[2] << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file.\n";
    }

    // Setup for potential energy profile (no plotting here)
    double Z_H_min = 2.8;
    double Z_H_max = 6.8;
    double step = 1e-2;
    std::vector<double> Z_H_vals;
    for (double z = Z_H_min; z <= Z_H_max; z += step) {
        Z_H_vals.push_back(z);
    }

    std::vector<double> V(Z_H_vals.size(), 0.0); // same as zeros(1, len)

    std::cout << "Setup complete. Atom count: " << Na << ", Z steps: " << Z_H_vals.size() << std::endl;

const double PI = 3.14159265358979323846;

double cosd(double angle) {
    return std::cos(angle * PI / 180.0);
}

double sind(double angle) {
    return std::sin(angle * PI / 180.0);
}

double max_val(const std::vector<double>& v) {
    return *std::max_element(v.begin(), v.end());
}

double min_val(const std::vector<double>& v) {
    return *std::min_element(v.begin(), v.end());
}

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

int main() {
    double pm = 2.47;
    int N = 15;
    double A = 15.2;
    double B = 24100;

    std::vector<double> a = {pm, 0.0, 0.0};
    std::vector<double> b = {pm * cosd(60), pm * sind(60), 0.0};
    std::vector<double> motif_1 = {0.0, 0.0, 0.0};
    std::vector<double> motif_2 = {pm / 2.0, pm / 2.0 * std::tan(PI / 6.0), 0.0};

    int Na = 2 * (N + 1) * (N + 1);
    std::vector<std::vector<double>> R_n(Na, std::vector<double>(3, 0.0));

    double d = pm / (2.0 * cosd(30));
    double D = std::sqrt(3.0 / 4.0) * d;

    int n = 0;
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (n < Na) {
                R_n[n][0] = i * a[0] + j * b[0] + motif_1[0];
                R_n[n][1] = i * a[1] + j * b[1] + motif_1[1];
                R_n[n][2] = i * a[2] + j * b[2] + motif_1[2];
                ++n;
            }
            if (n < Na) {
                R_n[n][0] = i * a[0] + j * b[0] + motif_2[0];
                R_n[n][1] = i * a[1] + j * b[1] + motif_2[1];
                R_n[n][2] = i * a[2] + j * b[2] + motif_2[2];
                ++n;
            }
        }
    }

    std::vector<double> X_s(Na), Y_s(Na);
    for (int i = 0; i < Na; ++i) {
        X_s[i] = R_n[i][0];
        Y_s[i] = R_n[i][1];
    }

    double X_H = max_val(X_s) / 2.0;
    double Y_H = max_val(Y_s) / 2.0;
    double X_T = X_H;
    double Y_T = Y_H + d;
    double X_B = X_H + D;
    double Y_B = Y_H;

    double Z_H_min = 2.8, Z_H_max = 6.8, step = 0.01;
    int steps = static_cast<int>((Z_H_max - Z_H_min) / step) + 1;
    std::vector<double> Z_H_vals = linspace(Z_H_min, Z_H_max, steps);

    std::vector<double> V(Z_H_vals.size(), 0.0);

    for (int i = 0; i < Na; ++i) {
        for (size_t j = 0; j < Z_H_vals.size(); ++j) {
            double Z_H = Z_H_vals[j];
            double R_H = std::sqrt(std::pow(X_H - X_s[i], 2) + std::pow(Y_H - Y_s[i], 2) + Z_H * Z_H);
            V[j] += B / std::pow(R_H, 12) - A / std::pow(R_H, 6);
        }
    }

    auto min_it = std::min_element(V.begin(), V.end());
    double Z_H_min_eq = Z_H_vals[std::distance(V.begin(), min_it)];
    std::cout << "Minimum energy at Z_H = " << Z_H_min_eq << " with value " << *min_it << std::endl;

    std::vector<double> V_B(Z_H_vals.size(), 0.0), V_T(Z_H_vals.size(), 0.0);
    for (int i = 0; i < Na; ++i) {
        for (size_t j = 0; j < Z_H_vals.size(); ++j) {
            double Z_H = Z_H_vals[j];
            double R_B = std::sqrt(std::pow(X_B - X_s[i], 2) + std::pow(Y_B - Y_s[i], 2) + Z_H * Z_H);
            double R_T = std::sqrt(std::pow(X_T - X_s[i], 2) + std::pow(Y_T - Y_s[i], 2) + Z_H * Z_H);
            V_B[j] += B / std::pow(R_B, 12) - A / std::pow(R_B, 6);
            V_T[j] += B / std::pow(R_T, 12) - A / std::pow(R_T, 6);
        }
    }

    double Z_H_ad = 3.38;
    std::vector<double> Xb = linspace(min_val(X_s), max_val(X_s), 500);
    std::vector<double> Yb = linspace(min_val(Y_s), max_val(Y_s), 500);

    int X_max = Xb.size(), Y_max = Yb.size();
    std::vector<std::vector<double>> V_c(Y_max, std::vector<double>(X_max, 0.0));

    for (int i = 0; i < Y_max; ++i) {
        for (int j = 0; j < X_max; ++j) {
            for (int k = 0; k < Na; ++k) {
                double R_C = std::sqrt(std::pow(Xb[j] - R_n[k][0], 2) + std::pow(Yb[i] - R_n[k][1], 2) + Z_H_ad * Z_H_ad);
                V_c[i][j] += B / std::pow(R_C, 12) - A / std::pow(R_C, 6);
            }
        }
    }

    std::cout << "Simulation completed." << std::endl;

  float bilinearInterp(const vector<float>& Xb, const vector<float>& Yb, const vector<vector<float>>& Vc, float x, float y) {
    int nx = Xb.size();
    int ny = Yb.size();

    int i = lower_bound(Xb.begin(), Xb.end(), x) - Xb.begin() - 1;
    int j = lower_bound(Yb.begin(), Yb.end(), y) - Yb.begin() - 1;

    if (i < 0 || j < 0 || i + 1 >= nx || j + 1 >= ny) return 0.0;

    float x1 = Xb[i], x2 = Xb[i + 1];
    float y1 = Yb[j], y2 = Yb[j + 1];
    float Q11 = Vc[j][i], Q21 = Vc[j][i + 1];
    float Q12 = Vc[j + 1][i], Q22 = Vc[j + 1][i + 1];

    float denom = (x2 - x1) * (y2 - y1);
    float val = 1.0 / denom * (
        Q11 * (x2 - x) * (y2 - y) +
        Q21 * (x - x1) * (y2 - y) +
        Q12 * (x2 - x) * (y - y1) +
        Q22 * (x - x1) * (y - y1)
    );

    return val;
}

int main() {
    // Assuming Xb, Yb, and Vc are already filled appropriately
    extern vector<float> Xb;
    extern vector<float> Yb;
    extern vector<vector<float>> Vc;

    vector<float> X_range, Y_range;
    float X_start = 29.5412, X_end = 32.0;
    float Y_start = 16.4, Y_end = 20.66;
    int Nx = 1000, Ny = 100;

    for (int i = 0; i < Nx; ++i) X_range.push_back(X_start + i * (X_end - X_start) / (Nx - 1));
    for (int i = 0; i < Ny; ++i) Y_range.push_back(Y_start + i * (Y_end - Y_start) / (Ny - 1));

    vector<float> V_c_X(Nx, 0.0f), V_c_Y(Ny, 0.0f);

    float Y_fixed = (Y_range.front() + Y_range.back()) / 2.0;
    float X_fixed = (X_range.front() + X_range.back()) / 2.0;

    for (int i = 0; i < Nx; ++i)
        V_c_X[i] = bilinearInterp(Xb, Yb, Vc, X_range[i], Y_fixed);

    for (int i = 0; i < Ny; ++i)
        V_c_Y[i] = bilinearInterp(Xb, Yb, Vc, X_fixed, Y_range[i]);

    // Optional: Write results to CSV for external plotting
    ofstream outX("Vc_X.csv");
    for (int i = 0; i < Nx; ++i)
        outX << X_range[i] << "," << V_c_X[i] << "\n";
    outX.close();

    ofstream outY("Vc_Y.csv");
    for (int i = 0; i < Ny; ++i)
        outY << Y_range[i] << "," << V_c_Y[i] << "\n";
    outY.close();

    cout << "Profiles exported to Vc_X.csv and Vc_Y.csv" << endl;
    return 0;
}
