#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include "include/runge_kutta4.hpp"

template<class T, size_t N>
static void save_data(const std::vector<T>& time, const std::vector<std::array<T, N>>& sol) {
    std::ofstream fout;
    for(size_t n = 0; n < N; ++n) {
        fout.open("X" + std::to_string(n) + ".csv");
        for(size_t i = 0; i < time.size(); ++i)
            fout << time[i] << ',' << sol[i][n] << std::endl;
        fout.close();
    }
}

int main() {
    static constexpr double A = 1.5, B = 3, V1 = 1, V2 = 1, l = 1;
    const auto [time, sol] =
        runge_kutta4(
            std::array<double, 2>{1, 1}, {0., 50}, 0.025,
            [](const double t, const std::array<double, 2>& x) noexcept {
                return std::array<double, 2>{
                    A - (B + 1) * x[0] + x[0]*x[0]*x[1],
                    B * x[0] - x[0]*x[0]*x[1]
                };
            }
        );

    save_data(time, sol);

    return 0;
}
