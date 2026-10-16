#ifndef trait_hpp__
#define trait_hpp__

#include "Population.hpp"
#include "dense_matrix.hpp"
#include "Utilites.hpp"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <filesystem>

namespace evogen
{
    class Trait
    {
    public:

        Trait();

#ifdef PYBIND

        Trait(Population &pop,
              pybind11::array_t<float> py_trmean,
              pybind11::array_t<float> py_qtl_prop_chrom,
              pybind11::array_t<float> py_corr_g,
              pybind11::array_t<float> py_varr_g,
              pybind11::array_t<float> py_corr_e,
              pybind11::array_t<float> py_varr_e,
              pybind11::array_t<float> py_envr,
              size_t dist_model,
              pybind11::array_t<float> py_dist_par);

        void set_trait(Population &pop,
                       pybind11::array_t<float> py_trmean,
                       pybind11::array_t<float> py_qtl_prop_chrom,
                       pybind11::array_t<float> py_corr_g,
                       pybind11::array_t<float> py_varr_g,
                       pybind11::array_t<float> py_corr_e,
                       pybind11::array_t<float> py_varr_e,
                       pybind11::array_t<float> py_envr,
                       size_t dist_model,
                       pybind11::array_t<float> py_dist_par);
        
        void reset_trait(Population &pop,
                       pybind11::array_t<float> py_trmean,
                       pybind11::array_t<float> py_qtl_prop_chrom,
                       pybind11::array_t<float> py_corr_g,
                       pybind11::array_t<float> py_varr_g,
                       pybind11::array_t<float> py_corr_e,
                       pybind11::array_t<float> py_varr_e,
                       pybind11::array_t<float> py_envr,
                       size_t dist_model,
                       pybind11::array_t<float> py_dist_par);

        void get_observations(Population &in_pop,
                              pybind11::array_t<float> py_env);

        void get_observations(Population &in_pop,
                              pybind11::array_t<float> py_env,
                              pybind11::array_t<float> py_t);
        
        void get_observations(Population &in_pop,
                              pybind11::array_t<float> py_env,                // reference to the base population for which making observations
                              pybind11::array_t<float> py_t, // container with observed trait values for each individual in in_pop
                              pybind11::array_t<int> py_g); // container with observed genotypes for each individual in in_pop

        void get_observations(Population &in_pop,
                              pybind11::array_t<float> py_env,
                              const std::string &out_t);
        
        void get_observations(Population &in_pop,
                              pybind11::array_t<float> py_env,
                              const std::string &out_t,  // output file name where observed trait values for each individual in in_pop will be writen
                              const std::string &out_g); // output file name where observed genotypes for each individual in in_pop will be writen
#endif
        Trait(Population &pop,
              std::vector<float> &trmean,
              std::vector<float> &qtl_prop_chrom,
              std::vector<std::vector<float>> &corr_g,
              std::vector<float> &varr_g,
              std::vector<std::vector<float>> &corr_e,
              std::vector<float> &varr_e,
              std::vector<float> &envr,
              size_t dist_model,
              std::vector<float> &dist_par);

        void set_trait(Population &pop,                          // reference to the base population
                       std::vector<float> &trmean,              // expected values of trait in base population, size=(n_trait,1)
                       std::vector<float> &qtl_prop_chrom,      // proportion of snps used as qtls for each chromosome, size=(n_chr,1)
                       std::vector<std::vector<float>> &corr_g, // correlations of genetic effects, size=(n_trait, n_trait)
                       std::vector<float> &varr_g,              // variances of genetic effects, size=(n_trait, 1)
                       std::vector<std::vector<float>> &corr_e, // correlations of environmental (residual) effects, size=(n_trait, n_trait)
                       std::vector<float> &varr_e,              // variances of environmental effects, size=(n_trait, 1)
                       std::vector<float> &envr,                // Control of an environment (accounting its changes), size=(n_trait,1)
                       size_t dist_model,                        // type of distribution used to sample dominance
                       std::vector<float> &dist_par);           // dominance distribution parameters, size=(2,1)

        void reset_trait(Population &pop,
                         std::vector<float> &trmean,
                         std::vector<float> &qtl_prop_chrom,
                         std::vector<std::vector<float>> &corr_g,
                         std::vector<float> &varr_g,
                         std::vector<std::vector<float>> &corr_e,
                         std::vector<float> &varr_e,
                         std::vector<float> &envr,
                         size_t dist_model,
                         std::vector<float> &dist_par);

        void get_observations(Population &in_pop,
                              std::vector<float> &env);

        void get_observations(Population &in_pop,
                              std::vector<float> &env,                // reference to the base population for which making observations
                              std::vector<std::vector<float>> &out_t, // container with observed trait values for each individual in in_pop
                              std::vector<std::vector<short>> &out_g); // container with observed genotypes for each individual in in_pop

        void get_observations(Population &in_pop,
                              std::vector<float> &env,
                              const std::string &out_t,  // output file name where observed trait values for each individual in in_pop will be writen
                              const std::string &out_g); // output file name where observed genotypes for each individual in in_pop will be writen

        void get_observations(Population &in_pop,
                              std::vector<float> &env,
                              std::vector<std::vector<float>> &out_t);

        void get_observations(Population &in_pop,
                              std::vector<float> &env,
                              const std::string &out_t);
        void clear();

        bool is_cleared();

        ~Trait();

        friend class Group;

    private:
        bool cleared;

        evolm::matrix<float> a;      // additive effects, size=(qtls,n_trait)
        evolm::matrix<float> e;      // environmental effects, size=(qtls,n_trait)
        evolm::matrix<float> k;      // dominance effects, size=(qtls,1)
        evolm::matrix<float> qtls;   // trait responsible qtls locations, size=(n_qtls,1)
        evolm::matrix<float> t_mean; // correction mean for the base population
        evolm::matrix<float> ta;     // genotype-determined trait
        evolm::matrix<float> te;     // environment-determined trait
        evolm::matrix<float> ae;     // additive genetic effect: ta - dom_effect

        std::vector<std::vector<unsigned long>> base_genome_structure;

        void sample_dom(size_t which_dist, std::vector<float> &dist_param);
        void sample_genes(std::vector<size_t> &n_qtl,
                          std::vector<std::vector<unsigned long>> &stable,
                          std::vector<float> &qtl_prop);
        void sample_effects(size_t n_trate);

        void calculate_trait(Population &in_pop, std::vector<float> &envr, size_t n_trate);
        void calculate_trait(Population &in_pop, std::vector<size_t> &ind_list, std::vector<float> &envr, size_t n_trate);

        float ploidy_effect(size_t n_ploidy, float degree);

        std::vector<int> get_locus_state(Population &in_pop,
                                         size_t individ,
                                         size_t qtl_position);

        evolm::matrix<float> var_diag(evolm::matrix<float> &arr);
        evolm::matrix<float> get_scaler(std::vector<float> &in_var, evolm::matrix<float> &in_arr);

        void realloc_traits(Population &in_pop, size_t n_trait);
        void realloc_traits(size_t pop_size, size_t n_trait);
        void calculate_correction_mean(std::vector<float> &in_mean);

#ifdef PYBIND
        std::vector<float> py_tovect(pybind11::array_t<float> py_vect);
        std::vector<std::vector<float>> py_tovect2d(pybind11::array_t<float> py_vect);
#endif

    protected:
    };

} // end of namespace evogen

#endif // trait_hpp__