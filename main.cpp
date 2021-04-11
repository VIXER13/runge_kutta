#include "include/runge_kutta.hpp"
#include <cmath>
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

    parameters<double> params;
    params.tau_init = 0.000001;
    params.time_interval = {t0, t0 + 4};
    const std::array<double, 2> init = {std::exp(std::sin(t0 * t0)), std::exp(std::cos(t0 * t0))};
    static constexpr auto system = [](const double t, const std::array<double, 2>& x) {
        return std::array<double, 2>{
             2 * t * x[0] * std::cos(t * t),
            -2 * t * x[1] * std::sin(t * t)
        };
    };

    auto [time, sol] = runge_kutta<double, runge_kutta4>(params, init, system);
    std::cout << "steps = " << time.size()-1 << std::endl;
    save_data("test", time, sol, 1000);

//    params.tau_init /= 2;
//    auto [time2, sol2] = runge_kutta<double, runge_kutta4>(params, std::array<double, 2>{1, 0}, system);
//    std::cout << "steps = " << time2.size()-1 << std::endl;
//
//    double err1 = 0, err2 = 0;
//    for(size_t i = 0; i < time.size()-1; ++i) {
//        err1 += std::abs(sol[i][0] - std::cos(time[i]));
//        err2 += std::abs(sol2[2*i][0] - std::cos(time[i]));
//    }
//    std::cout << "err1 = " << err1 << std::endl;
//    std::cout << "err2 = " << err2 << std::endl;
//    std::cout << "err1 / err2 = " << err1 / err2 << std::endl;

    return 0;
}
