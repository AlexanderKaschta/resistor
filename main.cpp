// Resistor - Tool to calculate resistor networks for a specific resistance out of E series resistors
// Copyright (C) 2024  Alexander Kaschta
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <chrono>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include <iomanip>


std::tuple<double, int> closest_simple(double searchValue, std::vector<double> entries) {
    // Find the closest item
    int closest_index = 0;
    double closest_value = entries[0];

    // I could start with 1 as 0th index is the default
    for (int i = 0; i < entries.size(); ++i) {
        if (std::abs(searchValue - entries[i]) < std::abs(searchValue - closest_value)) {
            closest_index = i;
            closest_value = entries[i];
        }
        // TODO: Check if abort condition (break if entries[i] > searchValue) will speed up the software
    }
    return std::make_tuple(closest_value, closest_index);
}

std::tuple<double, int, int> closest_series2(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    double closest_value = entries[0] + entries[0];

    int i_max = closest_index;
    if (searchValue > entries[closest_index] && closest_index < entries.size() - 1) {
        i_max = closest_index + 1;
    }

    for (int i = 0; i < i_max; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = entries[i] + entries[j];
            // If the new value is closer than the old value
            if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                // update the closest value
                closest_index_1 = i;
                closest_index_2 = j;
                closest_value = sum;
            }
            if (sum > searchValue) {
                break;
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2);
}

std::tuple<double, int, int> closest_parallel2(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    // if both entries are the same, the equation for parallel resistance reduces to a^2/(2a) or simply a/2, if a > 0.
    double closest_value = 0.5 * entries[0];

    for (int i = 0; i < entries.size(); ++i) {
        for (int j = 0; j < entries.size(); ++j) {
            double product = (entries[i] * entries[j]) / (entries[i] + entries[j]);
            // If the new value is closer than the old value
            if (std::abs(product - searchValue) < std::abs(closest_value - searchValue)) {
                // update the closest value
                closest_index_1 = i;
                closest_index_2 = j;
                closest_value = product;
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2);
}

std::tuple<double, int, int, int> closest_series3(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    double closest_value = 3.0 * entries[0];

    int i_max = closest_index;
    if (searchValue > entries[closest_index] && closest_index < entries.size() - 1) {
        i_max = closest_index + 1;
    }

    for (int i = 0; i < i_max; ++i) {
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k <= j; ++k) {
                double sum = entries[i] + entries[j] + entries[k];
                // If the new value is closer than the old value
                if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                    // update the closest value
                    closest_index_1 = i;
                    closest_index_2 = j;
                    closest_index_3 = k;
                    closest_value = sum;
                }
                if (sum > searchValue) {
                    break;
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3);
}

std::tuple<double, int, int, int>
closest_series3_parallel_in_series(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    double closest_value = 1.5 * entries[0];

    int i_max = closest_index;
    if (searchValue > entries[closest_index] && closest_index < entries.size() - 1) {
        i_max = closest_index + 1;
    }

    for (int i = 0; i <= i_max; ++i) {
        for (int j = 0; j < entries.size(); ++j) {
            for (int k = 0; k < entries.size(); ++k) {
                double sum = entries[i] + (entries[j] * entries[k]) / (entries[j] * entries[k]);
                // If the new value is closer than the old value
                if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                    // update the closest value
                    closest_index_1 = i;
                    closest_index_2 = j;
                    closest_index_3 = k;
                    closest_value = sum;
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3);
}

std::tuple<double, int, int, int>
closest_series3_series_in_parallel(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    double closest_value = ((entries[0] + entries[0]) * entries[0]) / (3.0 * entries[0]);


    for (int i = 0; i < entries.size(); ++i) {
        for (int j = 0; j < entries.size(); ++j) {
            for (int k = 0; k < entries.size(); ++k) {
                double sum = ((entries[i] + entries[j]) * entries[k]) / (entries[i] + entries[j] + entries[k]);
                // If the new value is closer than the old value
                if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                    // update the closest value
                    closest_index_1 = i;
                    closest_index_2 = j;
                    closest_index_3 = k;
                    closest_value = sum;
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3);
}

std::tuple<double, int, int, int>
closest_parallel3(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    double closest_value = (entries[0] * entries[0] * entries[0]) /
                           (entries[0] * entries[0] + entries[0] * entries[0] + entries[0] * entries[0]);

    for (int i = 0; i < entries.size(); ++i) {
        for (int j = 0; j < entries.size(); ++j) {
            for (int k = 0; k < entries.size(); ++k) {
                double sum = (entries[i] * entries[j] * entries[k]) /
                             (entries[j] * entries[k] + entries[i] * entries[k] + entries[i] * entries[j]);
                // If the new value is closer than the old value
                if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                    // update the closest value
                    closest_index_1 = i;
                    closest_index_2 = j;
                    closest_index_3 = k;
                    closest_value = sum;
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3);
}

std::tuple<double, int, int, int, int>
closest_series4(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    int closest_index_4 = 0;
    double closest_value = 4.0 * entries[0];

    int i_max = closest_index;
    if (searchValue > entries[closest_index] && closest_index < entries.size() - 1) {
        i_max = closest_index + 1;
    }

    for (int i = 0; i < i_max; ++i) {
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k <= j; ++k) {
                for (int l = 0; l <= k; ++l) {
                    double sum = entries[i] + entries[j] + entries[k] + entries[l];
                    // If the new value is closer than the old value
                    if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                        // update the closest value
                        closest_index_1 = i;
                        closest_index_2 = j;
                        closest_index_3 = k;
                        closest_index_4 = l;
                        closest_value = sum;
                    }
                    if (sum > searchValue) {
                        break;
                    }
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3, closest_index_4);
}

std::tuple<double, int, int, int, int>
closest_series4_parallel3_in_series(double searchValue, int closest_index, std::vector<double> entries) {
    int closest_index_1 = 0;
    int closest_index_2 = 0;
    int closest_index_3 = 0;
    int closest_index_4 = 0;
    double closest_value = entries[0] + (entries[0] * entries[0] * entries[0]) /
                                        (entries[0] * entries[0] + entries[0] * entries[0] + entries[0] * entries[0]);

    int i_max = closest_index;
    if (searchValue > entries[closest_index] && closest_index < entries.size() - 1) {
        i_max = closest_index + 1;
    }

    for (int i = 0; i < i_max; ++i) {
        for (int j = 0; j < entries.size(); ++j) {
            for (int k = 0; k < entries.size(); ++k) {
                for (int l = 0; l < entries.size(); ++l) {
                    double sum = entries[i] + (entries[j] * entries[k] * entries[l]) /
                                              (entries[k] * entries[l] + entries[j] * entries[l] +
                                               entries[j] * entries[k]);
                    // If the new value is closer than the old value
                    if (std::abs(sum - searchValue) < std::abs(closest_value - searchValue)) {
                        // update the closest value
                        closest_index_1 = i;
                        closest_index_2 = j;
                        closest_index_3 = k;
                        closest_index_4 = l;
                        closest_value = sum;
                    }
                }
            }
        }
    }

    return std::make_tuple(closest_value, closest_index_1, closest_index_2, closest_index_3, closest_index_4);
}

void search(double searchValue, std::vector<double> entries) {
    // Check for n = 1
    std::cout << "=====  n = 1  =====" << std::endl;
    const auto closest_1 = closest_simple(searchValue, entries);

    if (std::get<0>(closest_1) == searchValue) {
        std::cout << "Use a single resistor of " << std::get<0>(closest_1) << " Ohm!" << std::endl;
        // Abort as an optimal solution has been found
        return;
    }
    std::cout << "Closest index: " << std::get<1>(closest_1) << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_1) << std::endl;

    // Check for n = 2
    std::cout << "=====  n = 2  =====" << std::endl;

    const auto closest_2_series = closest_series2(searchValue, std::get<1>(closest_1), entries);

    if (std::get<0>(closest_2_series) == searchValue) {
        std::cout << "Use two resistors " << entries[std::get<1>(closest_2_series)] << " and "
                  << entries[std::get<2>(closest_2_series)] << " in series!" << std::endl;
        // Abort as an optimal solution has been found
        return;
    }

    const auto closest_2_parallel = closest_parallel2(searchValue, std::get<1>(closest_1), entries);
    if (std::get<0>(closest_2_parallel) == searchValue) {
        std::cout << "Use two resistors " << entries[std::get<1>(closest_2_parallel)] << " and "
                  << entries[std::get<2>(closest_2_parallel)] << " in parallel!" << std::endl;
        // Abort as an optimal solution has been found
        return;
    }

    std::cout << "Series:" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_2_series) << ", " << std::get<2>(closest_2_series)
              << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_2_series) << std::endl;

    std::cout << "Parallel:" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_2_parallel) << ", " << std::get<2>(closest_2_parallel)
              << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_2_parallel) << std::endl;

    // Check for n = 3
    std::cout << "=====  n = 3  =====" << std::endl;
    const auto closest_3_series = closest_series3(searchValue, std::get<1>(closest_1), entries);

    if (std::get<0>(closest_3_series) == searchValue) {
        std::cout << "Use three resistors " << entries[std::get<1>(closest_3_series)] << ", "
                  << entries[std::get<2>(closest_3_series)] << " and " << entries[std::get<3>(closest_3_series)]
                  << " in series!" << std::endl;
        return;
    }

    const auto closest_3_parallel_in_series = closest_series3_parallel_in_series(searchValue, std::get<1>(closest_1),
                                                                                 entries);
    if (std::get<0>(closest_3_parallel_in_series) == searchValue) {
        std::cout << "Some optimal solution found" << std::endl;
        return;
    }

    const auto closest_3_series_in_parallel = closest_series3_series_in_parallel(searchValue, std::get<1>(closest_1),
                                                                                 entries);

    if (std::get<0>(closest_3_series_in_parallel) == searchValue) {
        std::cout << "Some optimal solution found" << std::endl;
        return;
    }

    const auto closest_3_parallel = closest_parallel3(searchValue, std::get<1>(closest_1), entries);

    if (std::get<0>(closest_3_parallel) == searchValue) {
        std::cout << "Some optimal solution found" << std::endl;
        return;
    }

    std::cout << "Series:" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_3_series) << ", " << std::get<2>(closest_3_series)
              << ", " << std::get<3>(closest_3_series) << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_3_series) << std::endl;

    std::cout << "Parallel in series" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_3_parallel_in_series) << ", "
              << std::get<2>(closest_3_parallel_in_series)
              << ", " << std::get<3>(closest_3_parallel_in_series) << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_3_parallel_in_series) << std::endl;

    std::cout << "Series in parallel" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_3_series_in_parallel) << ", "
              << std::get<2>(closest_3_series_in_parallel)
              << ", " << std::get<3>(closest_3_series_in_parallel) << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_3_series_in_parallel) << std::endl;

    std::cout << "Parallel" << std::endl;
    std::cout << "Closest index: " << std::get<1>(closest_3_parallel) << ", " << std::get<2>(closest_3_parallel)
              << ", " << std::get<3>(closest_3_parallel) << std::endl;
    std::cout << "Closest value: " << std::setprecision(10) << std::get<0>(closest_3_parallel) << std::endl;

    // Check for n = 4
    std::cout << "=====  n = 4  =====" << std::endl;
    const auto closest_4_series = closest_series4(searchValue, std::get<1>(closest_1), entries);

    if (std::get<0>(closest_4_series) == searchValue) {
        std::cout << "Use four resistors " << entries[std::get<1>(closest_4_series)] << ", "
                  << entries[std::get<2>(closest_4_series)] << ", " << entries[std::get<3>(closest_4_series)] <<
                  " and " << entries[std::get<4>(closest_4_series)] << " in series!" << std::endl;
        return;
    }
}


int main() {
    double searchValue = 63059.0;

    std::cout << "Value to search: " << searchValue << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();

    // E3 series
    std::vector<double> series3{1.0, 2.2, 4.7};

    // E6 series
    std::vector<double> series6{1.0, 1.5, 2.2, 3.3, 4.7, 6.8};

    // E12 series
    std::vector<double> series12{1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2};

    // E24 series
    std::vector<double> series{1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7,
                               5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1};

    // E48 series
    std::vector<double> series48{1.00, 1.05, 1.10, 1.15, 1.21, 1.27, 1.33, 1.40, 1.47, 1.54, 1.62, 1.69, 1.78, 1.87,
                                 1.96, 2.05, 2.15, 2.26, 2.37, 2.49, 2.61, 2.74, 2.87, 3.01, 3.16, 3.32, 3.48, 3.65,
                                 3.83, 4.02, 4.22, 4.42, 4.64, 4.87, 5.11, 5.36, 5.62, 5.90, 6.19, 6.49, 6.81, 7.15,
                                 7.50, 7.87, 8.25, 8.66, 9.09, 9.53};

    // E96 series
    std::vector<double> series96{1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30, 1.33, 1.37,
                                 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87, 1.91,
                                 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32, 2.37, 2.43, 2.49, 2.55, 2.61, 2.67,
                                 2.74, 2.80, 2.87, 2.94, 3.01, 3.09, 3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74,
                                 3.83, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23,
                                 5.36, 5.49, 5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32,
                                 7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76};

    // E192 series
    std::vector<double> series192{1.00, 1.01, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.10, 1.11, 1.13, 1.14, 1.15, 1.17,
                                  1.18, 1.20, 1.21, 1.23, 1.24, 1.26, 1.27, 1.29, 1.30, 1.32, 1.33, 1.35, 1.37, 1.38,
                                  1.40, 1.42, 1.43, 1.45, 1.47, 1.49, 1.50, 1.52, 1.54, 1.56, 1.58, 1.60, 1.62, 1.64,
                                  1.65, 1.67, 1.69, 1.72, 1.74, 1.76, 1.78, 1.80, 1.82, 1.84, 1.87, 1.89, 1.91, 1.93,
                                  1.96, 1.98, 2.00, 2.03, 2.05, 2.08, 2.10, 2.13, 2.15, 2.18, 2.21, 2.23, 2.26, 2.29,
                                  2.32, 2.34, 2.37, 2.40, 2.43, 2.46, 2.49, 2.52, 2.55, 2.58, 2.61, 2.64, 2.67, 2.71,
                                  2.74, 2.77, 2.80, 2.84, 2.87, 2.91, 2.94, 2.98, 3.01, 3.05, 3.09, 3.12, 3.16, 3.20,
                                  3.24, 3.28, 3.32, 3.36, 3.40, 3.44, 3.48, 3.52, 3.57, 2.61, 3.65, 3.70, 3.74, 3.79,
                                  3.83, 3.88, 3.92, 3.97, 4.02, 4.07, 4.12, 4.17, 4.22, 4.27, 4.32, 4.37, 4.42, 4.48,
                                  4.53, 4.59, 4.64, 4.70, 4.75, 4.81, 4.87, 4.93, 4.99, 5.05, 5.11, 5.17, 5.23, 5.30,
                                  5.36, 5.42, 5.49, 5.56, 5.62, 5.69, 5.76, 5.83, 5.90, 5.97, 6.04, 6.12, 6.19, 6.26,
                                  6.34, 6.42, 6.49, 6.57, 6.65, 6.73, 6.81, 6.90, 6.98, 7.06, 7.15, 7.23, 7.32, 7.41,
                                  7.50, 7.59, 7.68, 7.77, 7.87, 7.79, 8.06, 8.16, 8.25, 8.35, 8.45, 8.56, 8.66, 8.76,
                                  8.87, 8.98, 9.09, 9.20, 9.31, 9.42, 9.53, 9.65, 9.76, 9.88};


    std::cout << "Number of entries in series: " << series.size() << std::endl;

    std::vector<double> power{10, 100, 1000, 10000, 100000};

    std::vector<double> entries;
    entries.reserve(power.size() * series.size());

    for (double pow: power) {
        for (double value: series) {
            entries.push_back(pow * value);
        }
    }

    std::cout << "Number of possible resistor values: " << entries.size() << std::endl;

    search(searchValue, entries);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "=====  TIME  =====" << std::endl;
    std::cout << "Duration in (ms): " << ms_double.count() << std::endl;

    return 0;
}
