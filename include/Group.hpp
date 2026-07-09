#ifndef group_hpp__
#define group_hpp__

#include "Population.hpp"
#include "Trait.hpp"
#include <vector>
#include <functional>
#include <algorithm>

namespace evogen
{
    class Group
    {
    public:
        Group();
        ~Group();

        size_t size();
        size_t size_at(size_t at);

        void clear();
        void remove();
        void move(Population &pop);
        void add(Population &pop);          // add entire population to the group
        void add(Population &pop,
                 size_t which_one);
        void add(Group &grp);
        void mate();                          // default: sexual reproduction, 1 - offspring with 1.0 success rate
        void mate(bool sexual_reproduction,   // true - sexual, false - asexsual
                  int max_offspring,          // max number of offspring from the same mating event
                  float success_rate,         // probability of getting max_offspring value (binomial distribution)
                  double mutation_frequency,  // mutation frequency - probability of mutation of a single snp during meyosis
                  size_t num_crossovers );
        void mate(bool sexual_reproduction,   // this is just to allow the float type for max_offspring variable (if someone will pass it in Python)
                  float max_offspring,
                  float success_rate,
                  double mutation_frequency,  // mutation frequency - probability of mutation of a single snp during meyosis
                  size_t num_crossovers );
        void regroup_newborn(Group &grp);

        void aging(int delta_t);
        void set_not_genotyped();
        void set_missing_observations();
        void kill();

        void select_basic( Group &selected, Group &notselected, const std::string category, const std::string mode, int num, float val = -1 );
        void select( const std::string category, const std::string mode, int num, float val = -1 );
        void select_into_group( Group &grp, const std::string category, const std::string mode, int num, float val = -1 );

        void get_genotypes(const std::string &file_out);
        void get_haplotypes(const std::string &file_out);
        void get_ancestry(const std::string &file_out);
        void get_pedigree(const std::string &file_out);
        void get_data(const std::string &file_out);
        void estimate_bv(const std::string &use_model, float var_g, float var_e);

        void get_genotypes_atlocus(size_t at_locus);

#ifdef PYBIND
        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              const std::string &out_trvalues);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              const std::string &out_trvalues,
                              const std::string &out_genotypes);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              pybind11::array_t<float> py_trvalues);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              pybind11::array_t<float> py_trvalues,
                              pybind11::array_t<int> py_genotypes);
#else
        void make_observation(Trait &trt,
                              std::vector<float> &env);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              const std::string &out_trvalues);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              const std::string &out_trvalues,
                              const std::string &out_genotypes);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              std::vector<std::vector<float>> &out_trvalues);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              std::vector<std::vector<float>> &out_trvalues,
                              std::vector<std::vector<short>> &out_genotypes);
#endif

    private:
        std::vector<std::reference_wrapper<Population>> pop_list;
        std::vector<std::vector<size_t>> individuals_list; // individuals_list[i][j] points to the position (index) in the active_individuals list

        std::vector<std::reference_wrapper<Population>> pop_list_newborn;
        std::vector<std::vector<size_t>> individuals_list_newborn;

        std::vector<size_t> get_individuals(size_t which_pop);
        void add_newborn(Population &pop, size_t which_one);

        // Comparison function to sort the vector elements by first element of tuples in descending order
        static bool sortdesc(const selection_candidate &a, const selection_candidate &b) { return (std::get<0>(a) > std::get<0>(b)); }

        std::vector<std::string> created_files; // the list of created files that needs to be deleted at the end (in destructor)
    };

} // end of namespace evogen

#endif // group_hpp__