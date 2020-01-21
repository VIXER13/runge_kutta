#ifndef DORMAND_PRINCE_HPP
#define DORMAND_PRINCE_HPP

#include <tuple>
#include <numeric>
#include <functional>

#include "array_extended.hpp"
#include "rows_different_sizes.hpp"

template<class Type>
struct parameters
{
    static_assert(std::is_floating_point<Type>::value, "The Type must be floating-point.");
    Type tol = 1e-9,            // Относительная точность вычислений
         factmax = 2.,          // Максимальный увеличивающий множитель
         factmin = 0.5,         // Минимальный увеличивающий множитель
         fact = 1.;             // Гарантийный множитель
    uint64_t min_segment = 100, // Минимальное число шагов
             max_segment = 500; // Максимальное число шагов
};

template<class Type>
rows_different_sizes<Type, uint8_t> init_butcher_table()
{
    rows_different_sizes<Type, uint8_t> table;
    table.reserve(22);
    table.push_back(std::array<Type, 1>({          0.                                                                 }));
    table.push_back(std::array<Type, 1>({    1./   5.                                                                 }));
    table.push_back(std::array<Type, 2>({    3./  40.,      9./  40.                                                  }));
    table.push_back(std::array<Type, 3>({   44./  45.,    -56./  15.,    32./   9.                                    }));
    table.push_back(std::array<Type, 4>({19372./6561., -25360./2187., 64448./6561., -212./729.                        }));
    table.push_back(std::array<Type, 5>({ 9017./3168.,   -355./  33., 46732./5247.,   49./176., -5103./18656.         }));
    table.push_back(std::array<Type, 6>({   35./ 384.,            0.,   500./1113.,  125./192., -2187./ 6784., 11./84.}));
    return table;
}

template<class Type, size_t Size>
std::tuple<std::vector<Type>, std::vector<std::array<Type, Size>>>
    dormand_prince(const std::array<Type, Size> &init, const std::array<Type, 2> time_interval, Type tau,
                   const std::function<std::array<Type, Size>(Type, std::array<Type, Size>)> &system, const parameters<Type> params = parameters<Type>())
{
    static const rows_different_sizes<Type, uint8_t> butcher_table_coeff  = init_butcher_table<double>();
    static const std::array<Type, 7>                 butcher_table_tau    = {          0., 1./5.,    3./   10.,   4./  5.,      8./     9.,         1.,     1.},
                                                     butcher_table_y_pred = {  35./  384.,    0.,  500./ 1113., 125./192.,  -2187./  6784.,  11./  84.,     0.},
                                                     butcher_table_y_corr = {5179./57600.,    0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.};

    Type tau_min = (time_interval[1]-time_interval[0]) / params.max_segment,
         tau_max = (time_interval[1]-time_interval[0]) / params.min_segment;

    std::cout << tau_min << " " << tau_max << std::endl;

    std::vector<Type> time = {time_interval[0]};
    std::vector<std::array<Type, Size>> sol = {init};

    std::array<Type, Size> y_pred = {}, y_corr = {}, K_pred = {}, K_corr;
    std::array<std::array<Type, Size>, 7> k = {};

    double err = 0.;
    while(time.back() < time_interval[1])
    {
        while(true)
        {
            for(size_t i = 0; i < 7; ++i)
            {
                K_pred = {};
                for(size_t j = 0; j < i; ++j)
                    K_pred += butcher_table_coeff(i, j) * k[j];
                k[i] = system(time.back() + butcher_table_tau[i]*tau, sol.back() + tau*K_pred);
            }

            K_pred = K_corr = {};
            for(size_t i = 0; i < 7; ++i)
            {
                K_pred += butcher_table_y_pred[i] * k[i];
                K_corr += butcher_table_y_corr[i] * k[i];
            }

            y_pred = sol.back() + tau*K_pred;
            y_corr = sol.back() + tau*K_corr;

            K_pred = y_pred - y_corr;

            // Обойдёмся простейшей нормой: суммой модулей
            err = std::accumulate(K_pred.cbegin(), K_pred.cend(), 0., [](const Type sum, const Type x) { return sum + fabs(x); }) /
                  std::accumulate(y_pred.cbegin(), y_pred.cend(), 0., [](const Type sum, const Type x) { return sum + fabs(x); }) / 31.;

            if(err < params.tol || tau < tau_min)
                break;
            else
                tau /= 2;
        }

        time.push_back(time.back()+tau);
        sol.push_back(y_pred);

        tau *= std::min(params.factmax, std::max(params.factmin, pow(params.tol/err, 1./6.)));
        tau = std::min(tau, tau_max);
    }

    return {time, sol};
}

#endif