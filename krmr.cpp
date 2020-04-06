#include <iostream>
#include <string>
#include <getopt.h>
#include <cstdint>
#include <functional>
#include "pliib.hpp"
//#include "kseq.hpp"
#include "spp.h"

struct hash_result{
    std::size_t num_hashes;
    std::uint64_t* hashes;
};

inline void kmer_to_integer(const char*& kmer, std::size_t kmer_len, std::int64_t& kmer_int){

}

inline void hash_sequence(const char*& sequence, std::size_t seq_len, std::uint8_t kmer_length, hash_result& hashes, auto& hash_func){
    
}


int main_filter(int argc, char** argv){
    std::vector<std::string> read_files;
    std::vector<std::string> ref_files;

    if (argc <= 2){
        std::cerr << "Please provide a GFA file." << std::endl;
        return -1;
    }

    optind = 2;
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"version", no_argument, 0, 'v'},

            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hpaAnel", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'f':
                read_files.push_back(optarg);
                break;
            case '?':
            case 'h':
                // nodes, edges, all stats, edges, paths
                exit(0);
            default:
                abort();
        }
    }

    /**
     * Read in an hash all the reference files.
     */

    return 0;

}
int main(int argc, char** argv){
    if (argc < 2){
        std::cerr << "No subcommand provided. Please provide a subcommand." << std::endl;
        return 1;
    }
     
    if (strcmp(argv[1], "filter") == 0){
        return main_filter(argc, argv);
    }
    else{
        std::cerr << "Invalid subcommand [" << argv[1] << "]. Please choose a valid subcommand (filter)" << std::endl;
    }
    return 0;
}
