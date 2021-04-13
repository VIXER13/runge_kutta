#include "include/runge_kutta.hpp"
#include <string_view>
#include <iostream>
#include <fstream>
#include <omp.h>

namespace {

template<class T, size_t N>
void save_data(const std::string& path,
               const std::vector<T>& time,
               const std::vector<std::array<T, N>>& sol,
               const std::vector<T>& loc_err,
               const size_t save_freq = 1) {
    std::ofstream fout;
    for(size_t i = 0; i < N; ++i) {
        fout.open(path + std::to_string(i) + ".csv");
        for(size_t j = 0; j < time.size(); ++j)
            if (j % save_freq == 0)
                fout << time[j] << ',' << sol[j][i] << '\n';
        fout.close();
    }
    fout.open(path + "_loc_err" + ".csv");
    for(size_t j = 0; j < time.size(); ++j)
        if (j % save_freq == 0)
            fout << time[j] << ',' << loc_err[j] << '\n';
    fout.close();
    fout.open(path + "_step" + ".csv");
    for(size_t j = 0; j < time.size()-1; ++j)
        if (j % save_freq == 0)
            fout << time[j] << ',' << time[j+1] - time[j] << '\n';
    fout.close();
}

}

int main() {
    static constexpr uintmax_t Nvar = 5;
    static constexpr double t0 = 0.1 * Nvar;
    static constexpr std::array<double, 2> init = {std::exp(std::sin(t0 * t0)),
                                                   std::exp(std::cos(t0 * t0))};
    static constexpr auto system = [](const double t, const std::array<double, 2>& x) {
        return std::array<double, 2>{
             2 * t * x[0] * std::log(std::max(x[1], 0.001)),
            -2 * t * x[1] * std::log(std::max(x[0], 0.001))
        };
    };
    ode::runge_kutta_parameters<double> parameters;
    parameters.time_interval = {t0, t0 + 4};
    parameters.tau_init = 0.1;
    parameters.tol = 1e-6;
    parameters.fact = 1;
    parameters.tol_small = parameters.tol * parameters.tol;
    double calc_time = omp_get_wtime();
    auto [time, sol, loc_err] = ode::runge_kutta<double, ode::runge_kutta4>(parameters, init, system);
    std::cout << "time = " << omp_get_wtime() - calc_time << std::endl;
    std::cout << "steps = " << time.size()-1 << std::endl;
    //save_data("eps3dp", time, sol, loc_err, 1);
    return EXIT_SUCCESS;
}
