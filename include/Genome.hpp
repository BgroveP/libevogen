#ifndef genome_hpp__
#define genome_hpp__

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <tuple>
#include "Utilites.hpp"
#include "Types_aliases.hpp"

namespace evogen
{
    class Genome
    {
    public:
        Genome();
        ~Genome();

        void set_genome(std::vector<std::vector<bool>> &snp,
                        std::vector<std::vector<unsigned long>> &gstructure,
                        std::vector<float> &gen_dist,
                        popid_t pop_origin_id);                          // uses prepared SNP variants either from file or from reproduction gamete; initializes gen_ancestry
        void set_genome(std::vector<std::vector<bool>> &snp,
                        std::vector<std::vector<ancestry_segment>> &ancestry,
                        std::vector<std::vector<unsigned long>> &gstructure,
                        std::vector<float> &gen_dist);               // uses prepared SNP variants either from file or from reproduction gamete; does not initialize gen_ancestry
        void set_genome(std::vector<std::vector<unsigned long>> &gstructure,
                        float ref_allele_probability,
                        short nploidy,
                        std::vector<float> &gen_dist,
                        popid_t pop_origin_id);                          // simulate SNP variants; initializes gen_ancestry
        
        void get_genome( std::vector<std::vector<bool>> &haplotypes,
                         std::vector<std::vector<unsigned long>> &gstructure ); /* returns genome information: variants and structure */
        void get_ancestry(std::vector<std::vector<ancestry_segment>> &ancestry_haplotypes);
        short get_ploidy();
        void get_reproduction_gamete(std::vector<std::vector<bool>> &out_gamete,
                                     std::vector< std::vector<ancestry_segment>> &out_ancestry,
                                     size_t cross_per_chr,
                                     double mut_freq);                   /* NOTE! we allow possibility for polyploid gametes */
        std::vector<short> get_genome();
        std::vector<short> get_genome_at(size_t locus);
        std::vector<std::vector<genlen_t>> get_snp_table();
        std::vector<std::vector<genlen_t>> get_genome_structure();
        std::vector<float> get_gen_distance();
        void clear();

        void operator=(const Genome& g)
        {
            markers = g.markers;
            structure = g.structure;
            snp_table = g.snp_table;
            hotspots = g.hotspots;
            std_hotspots = g.std_hotspots;
            gen_distance = g.gen_distance;
            gen_ancestry = g.gen_ancestry;
        }

        // For tests
        void show_genome( int max_variants );
        void show_recombination();


    private:
    
        std::vector<std::vector<bool>> markers;   // SNP (markers) variants for entire genome, the structure descriped in 'structure'; maternal comes first, than paternal
        std::vector<std::vector<genlen_t>> structure; // length of each chromosome (bp) and markers density (constant, bp)
        std::vector<std::vector<genlen_t>> snp_table; // for each chromosome indicates the first and the last snp in consecutive order
        std::vector<std::vector<size_t>> hotspots; // locations of recombination hotspots for each chromosome
        std::vector<float> std_hotspots; // allowed standart deviation of recombination events around specific hotspot
        std::vector<float> gen_distance; // genetic distance of each marker for all chromosomes (cM)
        std::vector<std::vector<ancestry_segment>> gen_ancestry; // vector of genomic segments (def. by tuples<gen.pos, pop_id>) marking their populational origin; has the same ploidy structure as markers vector; maternal comes first, than paternal

        //const float freq_hotspot = 0.05f; // frequency of recombination hotspots; n_hotspot = (freq_hotspot/mark_density) * n_snp_chr, number of equally distributed hotspots along each chromosome
        //const float step_devider = 5.0f;
        size_t n_hotspot = 10;

        short asign_snp_variant( float ref_allele_probability );

        // part of reproduction process
        void mutation( std::vector<std::vector<bool>> &in_gamete, double mut_freq );
        void recombination( std::vector<std::vector<bool>> &out_gamete, std::vector< std::vector<ancestry_segment>> &out_ancestry, size_t n_crosses );
        void crossover( size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide );
        void get_chromatide( size_t which_chr, size_t which_strand, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide );
        void get_ancestry_segment(size_t first_snp, size_t last_snp, size_t which_strand, size_t which_strands_set, std::vector<ancestry_segment> &ancestry_chromatide);
        
        void crossover_model_gamma(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide);
        void crossover_model_uniform(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide);
        void crossover_model_genmap(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide);
        double map_function_haldane(double distance);
        double map_function_kosambi(double distance);
        bool in_segment(std::vector<ancestry_segment> &which_segment, genlen_t which_item);
        void remove_redundant_segments(std::vector<ancestry_segment> &which_segment);
        void append_segment(std::vector<ancestry_segment> &which_segment, genlen_t which_item);
        int get_segment_pos(std::vector<ancestry_segment> &which_segment, genlen_t which_item);
        
        void def_snp_table();
        void def_hotspot_table();

    protected:
    };
} // end of namespace evo

#endif // genome_hpp__