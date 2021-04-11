#include "include/runge_kutta.hpp"
#include <string_view>
#include <iostream>
#include <fstream>

namespace {

template<class T, size_t N>
void save_data(const std::string& path, const std::vector<T>& time, const std::vector<std::array<T, N>>&sol,
               const size_t save_freq = 1) {
    std::ofstream fout;
    for(size_t i = 0; i < N; ++i) {
        fout.open(path + std::to_string(i) + ".csv");
        for(size_t j = 0; j < time.size(); ++j) {
            if (j % save_freq == 0)
                fout << time[j] << ',' << sol[j][i] << '\n';
        }
        fout.close();
    }
}

}

int main() {
    static constexpr uintmax_t Nvar = 5;
    static constexpr double t0 = 0.1 * Nvar;

    ode::runge_kutta_parameters<double> parameters;
    parameters.tau_init = 0.1;
    parameters.tol = 1e-9;
    parameters.time_interval = {t0, t0 + 4};
    const std::array<double, 2> init = {std::exp(std::sin(t0 * t0)), std::exp(std::cos(t0 * t0))};
    static constexpr auto system = [](const double t, const std::array<double, 2>& x) {
        return std::array<double, 2>{
             2 * t * x[0] * std::cos(t * t),
            -2 * t * x[1] * std::sin(t * t)
        };
    };

    auto [time, sol] = ode::runge_kutta<double, ode::dormand_prince>(parameters, init, system);
    std::cout << "steps = " << time.size()-1 << std::endl;
    save_data("test", time, sol, 1);

    return 0;
}
