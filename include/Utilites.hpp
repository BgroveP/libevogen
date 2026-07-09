#ifndef utilites_hpp__
#define utilites_hpp__

#include <iostream>
#include <random>
#include <chrono>
#include <filesystem>

namespace evogen
{
    class Utilites
    {
    public:
        Utilites();

        unsigned long get_randi(unsigned long range);

        template <typename T>
        T get_randi(T low_bound,
                    T upp_bound); /* returns a uniformly distributed random integer within the specified interval */

        template <typename T>
        std::vector<T> get_uni_rand(size_t n_values,
                                    T low_bound,
                                    T upp_bound,
                                    bool reset); /* returns a vector of uniformly distributed random integer within the specified interval */

        template <typename T>
        std::vector<T> get_runi_rand(size_t n_values,
                                    T low_bound,
                                    T upp_bound,
                                    bool reset); /* returns a vector of uniformly distributed random integer within the specified interval */

        template <typename T>
        std::vector<T> get_norm_rand(size_t n_values,
                                     T dist_mean,
                                     T dist_std,
                                     bool reset); /* returns a vector of normally distributed floating point random number */

        template <typename T>
        std::vector<T> get_gamma_rand(size_t n_values,
                                           T alpha,
                                           T beta,
                                           bool reset); /* returns a vector of gamma distributed random floating point values */

        std::vector<int> get_bin_rand(size_t n_values,
                                           int n,
                                           double p,
                                           bool reset); /* returns a vector of binary distributed random values */
                                           
        int bin_rand( int n, double p ); /* returns one binary distributed random value */
 
        void fremove(std::string file_name);

    private:
        // std::mt19937 generator;

    protected:
    };
} // end of namespace evogen
#endif // utilites_hpp__
