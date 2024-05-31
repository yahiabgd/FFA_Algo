#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include<chrono>
#include <functional>

#define M_PI       3.14159265358979323846   // pi



// Define the Firefly Algorithm (FFA) class
class FFA_V2 {
public:
    FFA_V2(int epoch, int pop_size, double gamma, double beta_base, double alpha, double alpha_damp, double delta, int exponent)
        : epoch_(epoch), pop_size_(pop_size), gamma_(gamma), beta_base_(beta_base), alpha_(alpha), alpha_damp_(alpha_damp),
        delta_(delta), exponent_(exponent) {}

    std::vector<double> solve(const std::function<double(const std::vector<double>&)>& obj_func,
        const std::vector<double>& lb, const std::vector<double>& ub) {
        // Initialize the population

        std::random_device rd;
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        std::vector<std::vector<double>> pop(pop_size_);
        for (auto& x : pop) {
            x.resize(lb.size());
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] = lb[i] + distribution(generator) * (ub[i] - lb[i]);
            }
        }
        std::vector<double> best_x(lb.size());
        double best_obj = std::numeric_limits<double>::max();

        // Main loop
        for (int e = 0; e < epoch_; ++e) {
            // Update the dynamic alpha
            alpha_ = alpha_damp_ * alpha_;

            // Evolve the population
            for (size_t i = 0; i < pop.size(); ++i) {
                std::vector<double> x = pop[i];
                std::vector<std::vector<double>> pop_child;
                for (size_t j = 0; j < pop.size(); ++j) {
                    // Calculate the radius and the attraction level
                    double dmax = std::sqrt(x.size());
                    double rij = 0.0;
                    for (size_t k = 0; k < x.size(); ++k) {
                        rij += std::pow(x[k] - pop[j][k], 2);
                    }
                    rij = std::sqrt(rij / (dmax * dmax));
                    double beta = beta_base_ * std::exp(-gamma_ * std::pow(rij, exponent_));

                    // Mutation vector
                    std::vector<double> mutation_vector(x.size());
                    for (size_t k = 0; k < x.size(); ++k) {
                        mutation_vector[k] = delta_ * (distribution(generator) - 0.5);
                    }

                    // Temp vector
                    std::vector<double> temp(x.size());
                    for (size_t k = 0; k < x.size(); ++k) {
                        temp[k] = (pop[j][k] - x[k]) * distribution(generator);
                    }

                    // New solution
                    std::vector<double> x_new(x.size());
                    for (size_t k = 0; k < x.size(); ++k) {
                        x_new[k] = x[k] + alpha_ * mutation_vector[k] + beta * temp[k];
                        // Boundary check
                        if (x_new[k] < lb[k]) x_new[k] = lb[k];
                        if (x_new[k] > ub[k]) x_new[k] = ub[k];
                    }

                    // Calculate the objective function value
                    double obj = obj_func(x_new);

                    // Update the population child
                    if (obj < obj_func(x)) {
                        pop_child.push_back(x_new);
                    }
                }

                // Add random solutions to the population child
                for (size_t j = pop_child.size(); j < pop_size_; ++j) {
                    std::vector<double> x_rand(x.size());
                    for (size_t k = 0; k < x.size(); ++k) {
                        x_rand[k] = lb[k] + distribution(generator) * (ub[k] - lb[k]);
                    }
                    pop_child.push_back(x_rand);
                }

                // Update the current best solution
                double obj_min = std::numeric_limits<double>::max();
                size_t idx_min = 0;
                for (size_t j = 0; j < pop_child.size(); ++j) {
                    double obj = obj_func(pop_child[j]);
                    if (obj < obj_min) {
                        obj_min = obj;
                        idx_min = j;
                    }
                }
                if (obj_min < obj_func(x)) {
                    pop[i] = pop_child[idx_min];
                }

                // Update the global best solution
                if (obj_min < best_obj) {
                    best_obj = obj_min;
                    best_x = pop_child[idx_min];
                }
            }

            // Print the current best solution
            std::cout << "Epoch " << e + 1 << ", best objective = " << best_obj<<std::endl ;//<< ", best solution = [";
            //for (size_t i = 0; i < best_x.size(); ++i) {
            //    std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
            //}
            //std::cout << std::endl;
        }
        for (size_t i = 0; i < best_x.size(); ++i) {
                std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
            }
            std::cout << std::endl;

        return best_x;
    }

private:
    int epoch_;
    int pop_size_;
    double gamma_;
    double beta_base_;
    double alpha_;
    double alpha_damp_;
    double delta_;
    int exponent_;
};

