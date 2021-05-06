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
    static constexpr double t0 = 0.1 * Nvar, tN = t0 + 4;
    static constexpr std::array init = {std::exp(std::sin(t0 * t0)),
                                        std::exp(std::cos(t0 * t0))};
    static constexpr auto system = [](const double t, const std::array<double, 2>& x) {
        return std::array{
             2 * t * x[0] * std::log(std::max(x[1], 1e-3)),
            -2 * t * x[1] * std::log(std::max(x[0], 1e-3))
        };
    };
    ode::runge_kutta_parameters<double> parameters;
    parameters.time_interval = {t0, tN};
    parameters.tol = 1e-4;
    for(const double tau : std::array{0.5, 0.1, 0.05})
        for(const double fact : std::array{0.9, 0.8})
            for(const double fact_max : std::array{1.5, 3.})
                for(const double fact_min : std::array{0.5, 0.2}) {
                    parameters.tau_init = tau;
                    parameters.fact = fact;
                    parameters.fact_max = fact_max;
                    parameters.fact_min = fact_min;
                    std::cout << "\\hline" << std::endl;
                    std::cout << tau << " & " << fact << " & " << fact_max << " & " << fact_min << " & ";
                    auto [time, sol, loc_err, def] = ode::runge_kutta<double, ode::runge_kutta4>(parameters, init, system);
                    std::cout << std::max(std::abs(sol.back()[0] - std::exp(std::sin(tN * tN))),
                                          std::abs(sol.back()[1] - std::exp(std::cos(tN * tN)))) << " & "
                              << def << " & " << time.size()-1 << " \\\\" << std::endl;
                }

    uintmax_t step = 0;
    parameters.fact = 0.9;
    parameters.fact_min = 0.2;
    parameters.tau_init = 0.2;
    for(const double tol : std::array{1e-3, 1e-4, 1e-6}) {
        parameters.tol = tol;
        parameters.tol_small = tol * tol;
        std::cout << "tol = " << tol << std::endl;
        {
            auto [time, sol, loc_err, def] = ode::runge_kutta<double, ode::runge_kutta4>(parameters, init, system);
            save_data("rk4"+std::to_string(step), time, sol, loc_err);
            std::cout << "rk4 steps = " << sol.size()-1 << " defect = " << def << " err = "
                      << std::max(std::abs(sol.back()[0] - std::exp(std::sin(tN * tN))),
                                  std::abs(sol.back()[1] - std::exp(std::cos(tN * tN)))) << std::endl;
        }
        {
            auto [time, sol, loc_err, def] = ode::runge_kutta<double, ode::dormand_prince>(parameters, init, system);
            save_data("dp"+std::to_string(step), time, sol, loc_err);
            std::cout << "dp steps = " << sol.size()-1 << " defect = " << def << " err = "
                      << std::max(std::abs(sol.back()[0] - std::exp(std::sin(tN * tN))),
                                  std::abs(sol.back()[1] - std::exp(std::cos(tN * tN)))) << std::endl;
        }
        ++step;
    }

    {
        parameters.fact = 0.9;
        parameters.fact_min = 0.2;
        parameters.tau_init = 0.2;
        parameters.tol = 1e-6;
        parameters.tau_init = 0.5;
        double t = omp_get_wtime();
        for(uintmax_t i = 0; i < 1000; ++i)
            auto [time, sol, loc_err, def] = ode::runge_kutta<double, ode::runge_kutta4>(parameters, init, system);
        std::cout << omp_get_wtime() - t << std::endl;

        t = omp_get_wtime();
        for(uintmax_t i = 0; i < 1000; ++i)
            auto [time, sol, loc_err, def] = ode::runge_kutta<double, ode::dormand_prince>(parameters, init, system);
        std::cout << omp_get_wtime() - t << std::endl;
    }

    return EXIT_SUCCESS;
}