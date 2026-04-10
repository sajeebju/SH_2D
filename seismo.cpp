#include "sh2dwave.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <string>

double get_env_double(const char* name, double default_value) {
    const char* val = std::getenv(name);
    return val ? std::stod(val) : default_value;
}

int get_env_int(const char* name, int default_value) {
    const char* val = std::getenv(name);
    return val ? std::stoi(val) : default_value;
}

void writeSeisToFile(const std::string& filename,
                     const std::pair<std::vector<double>, std::vector<double>>& data) {
    const auto& time = data.first;
    const auto& seis = data.second;

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::ios_base::failure("Failed to open file: " + filename);
    }

    for (size_t i = 0; i < time.size(); ++i) {
        outFile << time[i] << " " << seis[i] << "\n";
    }

    outFile.close();
    std::cout << "Data successfully written to " << filename << std::endl;
}

int main() {
    // Read parameters from environment
    double dx   = get_env_double("DX", 1.0);
    double dz   = get_env_double("DZ", 1.0);
    double dt   = get_env_double("DT", 0.001);
    double xmax = get_env_double("XMAX", 500.0);
    double zmax = get_env_double("ZMAX", 500.0);
    double tmax = get_env_double("TMAX", 0.502);

    double vs0  = get_env_double("VS0", 580.0);
    double rho0 = get_env_double("RHO0", 1000.0);

    double xsrc = get_env_double("XSRC", 250.0);
    double zsrc = get_env_double("ZSRC", 250.0);
    double xr   = get_env_double("XR", 330.0);
    double zr   = get_env_double("ZR", 330.0);

    double f0   = get_env_double("F0", 40.0);
    double t0   = get_env_double("T0", 4.0 / f0);

    int accuracy   = get_env_int("ACCURACY", 6);
    int use_absorb = get_env_int("USE_ABSORB", 1);
    int w          = get_env_int("W", 60);
    double a       = get_env_double("A", 0.0053);

    // Derived grid sizes
    int nx = static_cast<int>(xmax / dx);
    int nz = static_cast<int>(zmax / dz);
    int nt = static_cast<int>(tmax / dt);

    // Build homogeneous model
    std::vector<std::vector<double>> vs(nx, std::vector<double>(nz, vs0));
    std::vector<std::vector<double>> rho(nx, std::vector<double>(nz, rho0));

    // Convert physical coordinates to grid indices
    int isrc = static_cast<int>(xsrc / dx);
    int jsrc = static_cast<int>(zsrc / dz);
    int ir   = static_cast<int>(xr / dx);
    int jr   = static_cast<int>(zr / dz);

    auto result = SH_SEIS(nx, nz, nt, dx, dz, dt,
                          f0, t0,
                          isrc, jsrc, ir, jr,
                          vs, rho,
                          use_absorb, w, a,
                          accuracy);

    writeSeisToFile("seismogram.txt", result);
    std::cout << "Seismogram calculated successfully." << std::endl;

    return 0;
}