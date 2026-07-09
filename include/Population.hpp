#ifndef population_hpp__
#define population_hpp__

#include "Animal.hpp"
#include "Iointerface.hpp"
#include "Types_aliases.hpp"
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>
#include <map>
#include <bitset>
#include <limits>

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace evogen
{
    class Population
    {
    public:
        Population();
        ~Population();

        void set_population(size_t nindividuals,
                            const std::string &genotype_structure,
                            float ref_allele_probability,
                            int n_ploidy); /* set a 'wild' Population with simulated (or empty?) genome */

        void set_population(const std::string &haplotypes_fname,
                            const std::string &genotype_structure,
                            bool with_pedigree); /* sets Population with predefined genotypes which provided from file */

        void set_population(const std::string &plink_fname); // set population (with pedigree, sex, ids, family) from PLINK files
        
        size_t size();     // number of active individuals based on the active_individuals list
        size_t capacity(); // number of all (active + disabled) individuals, based on the individuals list
        void clear();
        void reshape();    // reduce capaciity to fit size (hence, memory)

        size_t get_popid();

        void aging(int delta_t);
        void set_birthday(int time);

        void id_at(size_t at, unsigned long id);
        unsigned long id_at(size_t at);

        void sire_at(size_t at, unsigned long id);
        unsigned long sire_at(size_t at);

        void dame_at(size_t at, unsigned long id);
        unsigned long dame_at(size_t at);

        void age_at(size_t at, int age);
        int age_at(size_t at);

        void birth_at(size_t at, int age);
        int birth_at(size_t at);

        void alive_at(size_t at, bool alive);
        bool alive_at(size_t at);
        
        void isgenotyped_at(size_t at, bool genotyped);
        bool isgenotyped_at(size_t at);

        void sex_at(size_t at, int sex);
        short sex_at(size_t at);

        short marker_phen_at(size_t at);

        size_t get_valid_pos(poplen_t in_pos);
        bool is_valid_pos(poplen_t in_pos);

#ifdef PYBIND
        void phen_at(size_t at, pybind11::array_t<float> phen); // do not make this public ?
        pybind11::array_t<float> phen_at(size_t at);

        void ebv_at(size_t at, pybind11::array_t<float> bv); // do not make this public ?
        pybind11::array_t<float> ebv_at(size_t at);

        void tbv_at(size_t at, pybind11::array_t<float> bv); // do not make this public ?
        pybind11::array_t<float> tbv_at(size_t at);

        void agv_at(size_t at, pybind11::array_t<float> bv); // do not make this public ?
        pybind11::array_t<float> agv_at(size_t at);

#endif
        void phen_at_cpp(size_t at, std::vector<float> &phen); // do not make this public ?
        std::vector<float> phen_at_cpp(size_t at);

        void ebv_at_cpp(size_t at, std::vector<float> &bv); // do not make this public ?
        std::vector<float> ebv_at_cpp(size_t at);

        void tbv_at_cpp(size_t at, std::vector<float> &bv); // do not make this public ?
        std::vector<float> tbv_at_cpp(size_t at);

        void agv_at_cpp(size_t at, std::vector<float> &bv); // do not make this public ?
        std::vector<float> agv_at_cpp(size_t at);

        // ------------- Not publik in Python ------------------------
        std::vector<short> get_genome_at(size_t which_genome, size_t locus);
        //void get_all_genotypes(const std::string &file_out);
        void get_all_genotypes(std::vector<std::vector<short>> &vect_out);
        void get_all_haplotypes(std::vector<std::vector<bool>> &vect_out,
                                std::vector<std::vector<unsigned long>> &out_snp_table,
                                std::vector<float> &out_gen_distance);
        
        void get_selected_haplotypes(std::vector<std::vector<bool>> &vect_out,
                                    std::vector<std::vector<unsigned long>> &out_snp_table,
                                    std::vector<float> &out_gen_distance,
                                    std::vector<poplen_t> &selected_individuals); // interface to access haplotypes from the Group class
        void get_selected_genotypes(const std::string &file_out, std::vector<poplen_t> &selected_individuals, bool append_mode); // interface to access genotypes from the Group class
        void get_selected_haplotypes(const std::string &file_out, std::vector<poplen_t> &selected_individuals, bool append_mode); // interface to access haplotypes from the Group class
        void get_selected_ancestry(const std::string &file_out, std::vector<poplen_t> &selected_individuals, bool append_mode); // interface to access ancestry from the Group class
        void get_selected_pedigree(const std::string &file_out, std::vector<poplen_t> &selected_individuals, bool append_mode); // interface to access pedigree from the Group class
        void get_selected_data(const std::string &file_out, std::vector<poplen_t> &selected_individuals, bool append_mode); // interface to access entire data (for LM) from the Group class

        void get_selected_genotypes_atlocus(std::vector<poplen_t> &selected_individuals, size_t at_locus); // interface to access genotypes from the Group class
        
        // ------------- Is publik in Python ------------------------
        void get_genotypes(const std::string &file_out);
        void get_haplotypes(const std::string &file_out);
        void get_ancestry(const std::string &file_out);
        void get_ld(const std::string &out_file,
                    bool full_info = false,
                    int which_chr = -1,
                    unsigned int snp_step = 1); // can be accessed only from Population class!
        void get_pedigree(const std::string &file_out);
        void get_data(const std::string &file_out);

        void get_genotypes_atlocus(size_t at_locus); // for every individual assign genotype value at a specified locus to an animal's property marker_phen

        // ------------- Required in Trait class ---------------------
        std::vector<std::vector<unsigned long>> get_genome_table();
        std::vector<std::vector<unsigned long>> get_genome_structure();
        size_t get_ploidy();
        // ------------------------------------------------------------

        friend class Group;

        // ------------------ TESTING --------------------------------
#ifdef UTEST
        void show_animals(size_t max_animals, size_t max_snps);

        void show_pop()
        {
            for (size_t i = 0; i < individuals.size(); i++)
            {
                std::cout<<"from individuals: i = "<<i
                         <<", id = "<<individuals[i].get_id()
                         <<", active = "<<individuals[i].get_active()
                         <<", sex = "<<individuals[i].get_sex()
                         <<", sire = "<<individuals[i].get_sire()
                         <<", dame = "<<individuals[i].get_dame()
                         <<", age = "<<individuals[i].get_age()
                         <<"\n";
            }
            std::cout<<"\n";
        }
#endif

    private:
        std::vector<Animal> individuals;
        std::vector<poplen_t> active_individuals;
        popid_t origin_id = 0; // pop id used to track genomic ancestry

        size_t last_added_index = 0; // track a last index in active_individuals during a last update of the vector
        size_t contin_individ_index = 0; // continiouslt tracking the number of added new individuals regardless of reshaping

        const std::string gdistance_name = "Genetic_Map(cM)";
        const std::string chr_name = "chr";

        void remove_at(size_t which_list_position);
        void add(Animal &in_a); // adds existent entiety from another population

        void add_at(Animal &in_a, size_t position);
        void resize(size_t n_elements);
        popid_t assign_origin_id();
    };

    inline bool operator==(const Population &rhs, const Population &lhs)
    {
        return std::addressof(rhs) == std::addressof(lhs);
    }

    inline bool operator!=(const Population &rhs, const Population &lhs)
    {
        return std::addressof(rhs) != std::addressof(lhs);
    }

} // end of namespace evogen

#endif // population_hpp__