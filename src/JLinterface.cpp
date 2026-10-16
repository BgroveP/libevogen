#include "jlcxx.hpp"
#include "Population.hpp"
#include "Group.hpp"
#include "Types_aliases.hpp"
#include "stl.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module &mod)
{
        // Adding types for nested vectors
        //    mod.add_type<std::vector<std::vector<float>>>("Vector{Vector{Float32}}")
        //       .constructor<>();
        //   mod.add_type<std::vector<std::vector<short>>>("Vector{Vector{Int16}}")
        //       .constructor<>();
        //   mod.add_type<std::vector<std::vector<unsigned long>>>("Vector{Vector{UInt64}}")
        //       .constructor<>();
        //   mod.add_type<std::vector<std::vector<bool>>>("Vector{Vector{Bool}}")
        //       .constructor<>();

        // Adding the Population Class
        mod.add_type<evogen::Population>("Population")
            .constructor<>()
            // --- Set population methods ---
            .method("set_population",
                    static_cast<void (evogen::Population::*)(
                        size_t,
                        const std::string &,
                        float,
                        int)>(&evogen::Population::set_population))
            .method("set_population",
                    static_cast<void (evogen::Population::*)(
                        const std::string &,
                        const std::string &,
                        bool)>(&evogen::Population::set_population))
            .method("set_population",
                    static_cast<void (evogen::Population::*)(
                        const std::string &)>(&evogen::Population::set_population))

            // --- Utility functions ---
            .method("size",
                    static_cast<size_t (evogen::Population::*)()>(&evogen::Population::size))
            .method("capacity",
                    static_cast<size_t (evogen::Population::*)()>(&evogen::Population::capacity))
            .method("clear",
                    static_cast<void (evogen::Population::*)()>(&evogen::Population::clear))
            .method("reshape",
                    static_cast<void (evogen::Population::*)()>(&evogen::Population::reshape))
            .method("get_popid",
                    static_cast<size_t (evogen::Population::*)()>(&evogen::Population::get_popid))

            // --- Age/Birthday ---
            .method("aging",
                    static_cast<void (evogen::Population::*)(int)>(&evogen::Population::aging))
            .method("set_birthday",
                    static_cast<void (evogen::Population::*)(int)>(&evogen::Population::set_birthday))

            // --- ID/Sire/Dame ---
            .method("id_at",
                    static_cast<void (evogen::Population::*)(size_t, unsigned long)>(&evogen::Population::id_at))
            .method("id_at",
                    static_cast<unsigned long (evogen::Population::*)(size_t)>(&evogen::Population::id_at))
            .method("sire_at",
                    static_cast<void (evogen::Population::*)(size_t, unsigned long)>(&evogen::Population::sire_at))
            .method("sire_at",
                    static_cast<unsigned long (evogen::Population::*)(size_t)>(&evogen::Population::sire_at))
            .method("dame_at",
                    static_cast<void (evogen::Population::*)(size_t, unsigned long)>(&evogen::Population::dame_at))
            .method("dame_at",
                    static_cast<unsigned long (evogen::Population::*)(size_t)>(&evogen::Population::dame_at))

            // --- Age/Birth/Alive/Genotyped ---
            .method("age_at",
                    static_cast<void (evogen::Population::*)(size_t, int)>(&evogen::Population::age_at))
            .method("age_at",
                    static_cast<int (evogen::Population::*)(size_t)>(&evogen::Population::age_at))
            .method("birth_at",
                    static_cast<void (evogen::Population::*)(size_t, int)>(&evogen::Population::birth_at))
            .method("birth_at",
                    static_cast<int (evogen::Population::*)(size_t)>(&evogen::Population::birth_at))
            .method("alive_at",
                    static_cast<void (evogen::Population::*)(size_t, bool)>(&evogen::Population::alive_at))
            .method("alive_at",
                    static_cast<bool (evogen::Population::*)(size_t)>(&evogen::Population::alive_at))
            .method("isgenotyped_at",
                    static_cast<void (evogen::Population::*)(size_t, bool)>(&evogen::Population::isgenotyped_at))
            .method("isgenotyped_at",
                    static_cast<bool (evogen::Population::*)(size_t)>(&evogen::Population::isgenotyped_at))

            // --- Sex/Marker Phenotype ---
            .method("sex_at",
                    static_cast<void (evogen::Population::*)(size_t, int)>(&evogen::Population::sex_at))
            .method("sex_at",
                    static_cast<short (evogen::Population::*)(size_t)>(&evogen::Population::sex_at))
            .method("marker_phen_at",
                    static_cast<short (evogen::Population::*)(size_t)>(&evogen::Population::marker_phen_at))

            // --- Genome/Genotype/Haplotype Accessors ---
            .method("get_genome_at",
                    static_cast<std::vector<short> (evogen::Population::*)(size_t, size_t)>(&evogen::Population::get_genome_at))
            .method("get_genome_table",
                    static_cast<std::vector<std::vector<unsigned long>> (evogen::Population::*)()>(&evogen::Population::get_genome_table))
            .method("get_genome_structure",
                    static_cast<std::vector<std::vector<unsigned long>> (evogen::Population::*)()>(&evogen::Population::get_genome_structure))
            .method("get_ploidy",
                    static_cast<size_t (evogen::Population::*)()>(&evogen::Population::get_ploidy))

            // --- Phenotype/EBV/TBV/AGV (C++ and Python) ---
            // C++ versions (std::vector<float>)
            .method("phen_at_cpp",
                    static_cast<void (evogen::Population::*)(size_t, std::vector<float> &)>(&evogen::Population::phen_at_cpp))
            .method("phen_at_cpp",
                    static_cast<std::vector<float> (evogen::Population::*)(size_t)>(&evogen::Population::phen_at_cpp))
            .method("ebv_at_cpp",
                    static_cast<void (evogen::Population::*)(size_t, std::vector<float> &)>(&evogen::Population::ebv_at_cpp))
            .method("ebv_at_cpp",
                    static_cast<std::vector<float> (evogen::Population::*)(size_t)>(&evogen::Population::ebv_at_cpp))
            .method("tbv_at_cpp",
                    static_cast<void (evogen::Population::*)(size_t, std::vector<float> &)>(&evogen::Population::tbv_at_cpp))
            .method("tbv_at_cpp",
                    static_cast<std::vector<float> (evogen::Population::*)(size_t)>(&evogen::Population::tbv_at_cpp))
            .method("agv_at_cpp",
                    static_cast<void (evogen::Population::*)(size_t, std::vector<float> &)>(&evogen::Population::agv_at_cpp))
            .method("agv_at_cpp",
                    static_cast<std::vector<float> (evogen::Population::*)(size_t)>(&evogen::Population::agv_at_cpp))

            // --- File Output Methods ---
            .method("get_genotypes",
                    static_cast<void (evogen::Population::*)(const std::string &)>(&evogen::Population::get_genotypes))
            .method("get_haplotypes",
                    static_cast<void (evogen::Population::*)(const std::string &)>(&evogen::Population::get_haplotypes))
            .method("get_ancestry",
                    static_cast<void (evogen::Population::*)(const std::string &)>(&evogen::Population::get_ancestry))
            .method("get_ld",
                    static_cast<void (evogen::Population::*)(const std::string &, bool, int, unsigned int)>(&evogen::Population::get_ld))
            .method("get_pedigree",
                    static_cast<void (evogen::Population::*)(const std::string &)>(&evogen::Population::get_pedigree))
            .method("get_data",
                    static_cast<void (evogen::Population::*)(const std::string &)>(&evogen::Population::get_data))
            .method("get_genotypes_atlocus",
                    static_cast<void (evogen::Population::*)(size_t)>(&evogen::Population::get_genotypes_atlocus))

            // --- Selected Data Accessors ---
            .method("get_selected_genotypes",
                    static_cast<void (evogen::Population::*)(const std::string &, std::vector<poplen_t> &, bool)>(&evogen::Population::get_selected_genotypes))
            .method("get_selected_haplotypes",
                    static_cast<void (evogen::Population::*)(const std::string &, std::vector<poplen_t> &, bool)>(&evogen::Population::get_selected_haplotypes))
            .method("get_selected_ancestry",
                    static_cast<void (evogen::Population::*)(const std::string &, std::vector<poplen_t> &, bool)>(&evogen::Population::get_selected_ancestry))
            .method("get_selected_pedigree",
                    static_cast<void (evogen::Population::*)(const std::string &, std::vector<poplen_t> &, bool)>(&evogen::Population::get_selected_pedigree))
            .method("get_selected_data",
                    static_cast<void (evogen::Population::*)(const std::string &, std::vector<poplen_t> &, bool)>(&evogen::Population::get_selected_data))
            .method("get_selected_genotypes_atlocus",
                    static_cast<void (evogen::Population::*)(std::vector<poplen_t> &, size_t)>(&evogen::Population::get_selected_genotypes_atlocus))

            // --- Vector Output Methods (for internal use) ---
            .method("get_all_genotypes",
                    static_cast<void (evogen::Population::*)(std::vector<std::vector<short>> &)>(&evogen::Population::get_all_genotypes))
            .method("get_all_haplotypes",
                    static_cast<void (evogen::Population::*)(
                        std::vector<std::vector<bool>> &,
                        std::vector<std::vector<unsigned long>> &,
                        std::vector<float> &)>(&evogen::Population::get_all_haplotypes))
            .method("get_selected_haplotypes",
                    static_cast<void (evogen::Population::*)(
                        std::vector<std::vector<bool>> &,
                        std::vector<std::vector<unsigned long>> &,
                        std::vector<float> &,
                        std::vector<poplen_t> &)>(&evogen::Population::get_selected_haplotypes))

            // --- Valid Position Check ---
            .method("get_valid_pos",
                    static_cast<size_t (evogen::Population::*)(poplen_t)>(&evogen::Population::get_valid_pos))
            .method("is_valid_pos",
                    static_cast<bool (evogen::Population::*)(poplen_t)>(&evogen::Population::is_valid_pos));

        // Register reference wrapper for Population
        mod.add_type<std::reference_wrapper<evogen::Population>>("Ref{Population}")
            .constructor<evogen::Population &>();

        // ======================
        // Trait Class Bindings
        // ======================
        mod.add_type<evogen::Trait>("Trait")
            .constructor<>() // Default constructor

            // Parameterized constructor
            .constructor<evogen::Population &,
                         std::vector<float> &,
                         std::vector<float> &,
                         std::vector<std::vector<float>> &,
                         std::vector<float> &,
                         std::vector<std::vector<float>> &,
                         std::vector<float> &,
                         std::vector<float> &,
                         size_t,
                         std::vector<float> &>()

            // Set and reset trait methods
            .method("set_trait",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &,
                                                        std::vector<float> &,
                                                        std::vector<float> &,
                                                        size_t,
                                                        std::vector<float> &)>(&evogen::Trait::set_trait))

            .method("reset_trait",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &,
                                                        std::vector<float> &,
                                                        std::vector<float> &,
                                                        size_t,
                                                        std::vector<float> &)>(&evogen::Trait::reset_trait))

            // Get observations methods
            .method("get_observations",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &)>(&evogen::Trait::get_observations))

            .method("get_observations",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &,
                                                        std::vector<std::vector<short>> &)>(&evogen::Trait::get_observations))

            .method("get_observations",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        const std::string &)>(&evogen::Trait::get_observations))

            .method("get_observations",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        const std::string &,
                                                        const std::string &)>(&evogen::Trait::get_observations))

            .method("get_observations",
                    static_cast<void (evogen::Trait::*)(evogen::Population &,
                                                        std::vector<float> &,
                                                        std::vector<std::vector<float>> &)>(&evogen::Trait::get_observations))

            // Clear method
            .method("clear",
                    static_cast<void (evogen::Trait::*)()>(&evogen::Trait::clear))

            // Is cleared method
            .method("is_cleared",
                    static_cast<bool (evogen::Trait::*)()>(&evogen::Trait::is_cleared));

        // ======================
        // Group Class Bindings
        // ======================
        mod.add_type<evogen::Group>("Group")
            .constructor<>()

            // Size methods
            .method("size",
                    static_cast<size_t (evogen::Group::*)()>(&evogen::Group::size))
            .method("size_at",
                    static_cast<size_t (evogen::Group::*)(size_t)>(&evogen::Group::size_at))

            // Clear and remove methods
            .method("clear",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::clear))
            .method("remove",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::remove))

            // Population management methods
            .method("move",
                    static_cast<void (evogen::Group::*)(evogen::Population &)>(&evogen::Group::move))
            .method("add",
                    static_cast<void (evogen::Group::*)(evogen::Population &)>(&evogen::Group::add))
            .method("add",
                    static_cast<void (evogen::Group::*)(evogen::Population &, size_t)>(&evogen::Group::add))
            .method("add",
                    static_cast<void (evogen::Group::*)(evogen::Group &)>(&evogen::Group::add))

            // Mate methods
            .method("mate",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::mate))
            .method("mate",
                    static_cast<void (evogen::Group::*)(bool, int, float, double, size_t)>(&evogen::Group::mate))
            .method("mate",
                    static_cast<void (evogen::Group::*)(bool, float, float, double, size_t)>(&evogen::Group::mate))

            // Regroup newborns
            .method("regroup_newborn",
                    static_cast<void (evogen::Group::*)(evogen::Group &)>(&evogen::Group::regroup_newborn))

            // Aging and status methods
            .method("aging",
                    static_cast<void (evogen::Group::*)(int)>(&evogen::Group::aging))
            .method("set_not_genotyped",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::set_not_genotyped))
            .method("set_missing_observations",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::set_missing_observations))
            .method("kill",
                    static_cast<void (evogen::Group::*)()>(&evogen::Group::kill))

            // Selection methods
            .method("select_basic",
                    static_cast<void (evogen::Group::*)(evogen::Group &, evogen::Group &, const std::string, const std::string, int, float)>(&evogen::Group::select_basic))
            .method("select",
                    static_cast<void (evogen::Group::*)(const std::string, const std::string, int, float)>(&evogen::Group::select))
            .method("select_into_group",
                    static_cast<void (evogen::Group::*)(evogen::Group &, const std::string, const std::string, int, float)>(&evogen::Group::select_into_group))

            // File output methods
            .method("get_genotypes",
                    static_cast<void (evogen::Group::*)(const std::string &)>(&evogen::Group::get_genotypes))
            .method("get_haplotypes",
                    static_cast<void (evogen::Group::*)(const std::string &)>(&evogen::Group::get_haplotypes))
            .method("get_ancestry",
                    static_cast<void (evogen::Group::*)(const std::string &)>(&evogen::Group::get_ancestry))
            .method("get_pedigree",
                    static_cast<void (evogen::Group::*)(const std::string &)>(&evogen::Group::get_pedigree))
            .method("get_data",
                    static_cast<void (evogen::Group::*)(const std::string &)>(&evogen::Group::get_data))

            // Estimate breeding values
            .method("estimate_bv",
                    static_cast<void (evogen::Group::*)(const std::string &, float, float)>(&evogen::Group::estimate_bv))

            // Genotypes at locus
            .method("get_genotypes_atlocus",
                    static_cast<void (evogen::Group::*)(size_t)>(&evogen::Group::get_genotypes_atlocus))

            // Make observation methods (C++ versions)
            .method("make_observation",
                    static_cast<void (evogen::Group::*)(evogen::Trait &, std::vector<float> &)>(&evogen::Group::make_observation))
            .method("make_observation",
                    static_cast<void (evogen::Group::*)(evogen::Trait &, std::vector<float> &, const std::string &)>(&evogen::Group::make_observation))
            .method("make_observation",
                    static_cast<void (evogen::Group::*)(evogen::Trait &, std::vector<float> &, const std::string &, const std::string &)>(&evogen::Group::make_observation))
            .method("make_observation",
                    static_cast<void (evogen::Group::*)(evogen::Trait &, std::vector<float> &, std::vector<std::vector<float>> &)>(&evogen::Group::make_observation))
            .method("make_observation",
                    static_cast<void (evogen::Group::*)(evogen::Trait &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)>(&evogen::Group::make_observation));

}
