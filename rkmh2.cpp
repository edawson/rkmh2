#include <iostream>
#include <string>
#include <getopt.h>
#include <sstream>
#include <cstdint>
#include <functional>
#include "pliib.hpp"
//#include "kseq.hpp"
#include "spp.h"
#include "mkmh.hpp"
#include "seqio.hpp"

//#define DEBUG_KRMR

struct hash_result{
    std::size_t num_hashes;
    mkmh::hash_t* hashes = nullptr;
    hash_result(){

    }
    hash_result(const std::size_t& hashcount){
        num_hashes = hashcount;
        hashes = new mkmh::hash_t [hashcount];
    }
    ~hash_result(){
        //num_hashes = 0;
        //delete [] hashes;
    }
    void clear(){
        num_hashes = 0;
        delete [] hashes;
    }
    void sort(){
        std::sort(this->hashes, this->hashes + this->num_hashes);
    }
    void sketch(std::size_t sketch_size = 1000){
        sort();
        std::size_t first_valid_index = 0;
        while (hashes[first_valid_index] == 0){
            ++first_valid_index;
        }
        std::size_t val_count = num_hashes - first_valid_index;
        std::size_t s = val_count < sketch_size ? val_count : sketch_size;
        mkmh::hash_t* sketch = new mkmh::hash_t[s];
        for (std::size_t i = first_valid_index; i < s; ++i){
            sketch[i] = hashes[i];
        }
        delete [] hashes;
        this->hashes = sketch;
        this->num_hashes = s;
    }
};

inline void hash_kmer(char* kmer,
        const std::size_t& kmer_len,
        mkmh::hash_t& ret,
        bool drop_amb_nucs = true){

    if (drop_amb_nucs){
        if (!pliib::canonical(kmer, kmer_len)){
            ret = 0;
            return;
        }
    }
   
    uint32_t rhash[4];
    uint32_t fhash[4];
    
    mkmh::mmh3::MurmurHash3_x64_128(kmer, kmer_len, 42, fhash);
    pliib::reverse_inplace(kmer, kmer_len);
    mkmh::mmh3::MurmurHash3_x64_128(kmer, kmer_len, 42, rhash);
    pliib::reverse_inplace(kmer, kmer_len);

    mkmh::hash_t tmp_fwd = static_cast<uint64_t>(fhash[0]) << 32 | fhash[1];
    mkmh::hash_t tmp_rev = static_cast<uint64_t>(rhash[0]) << 32 | rhash[1];
    ret = tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev;
}

inline void hash_sequence(char*& sequence,
        const std::size_t& seq_len,
        const int& kmer_length,
        hash_result& hashes,
        bool drop_amb_nucs = true,
        bool enforce_capital_nucs = true){
    
    if (enforce_capital_nucs){
        pliib::to_upper(sequence, seq_len);
    }
 
   
    std::size_t num_kmers = seq_len - kmer_length;
    hashes.num_hashes = num_kmers;
    hashes.hashes = new mkmh::hash_t [num_kmers];
    //#pragma omp for
    for (std::size_t i = 0; i < num_kmers; ++i){
        hash_kmer(sequence + i, kmer_length, hashes.hashes[i], drop_amb_nucs);
    }
}

struct read_t{
    char* seq = nullptr;
    char* name = nullptr;
    char* qual = nullptr;
    char* comment = nullptr;
    std::size_t seqlen = 0;
    read_t(){

    }
    ~read_t(){
        //delete [] seq;
        //delete [] name;
        //delete [] qual;
        //delete [] comment;
    }

    void purge_seq(){
        delete [] seq;
    }
    void clear(){
        delete [] seq;
        delete [] name;
        delete [] qual;
        delete [] comment;
    }

};

std::ostream& output_kseq(klibpp::KSeq& ks, std::ostream& os){
    if (ks.qual.empty()){
        os << ">" << ks.name << std::endl << ks.seq;
    }
    else{
        os << "@" << ks.name << std::endl <<
            ks.seq << std::endl <<
            "+" << ks.comment << std::endl <<
            ks.qual;
    }
    return os;
}

void usage(char** argv){
    std::cout << argv[1] << ": Filter reads using MinHash." << std::endl;
    std::cout << "Usage:" << argv[1] << " <subcommand> -r <ref FASTA/FASTQ> -f <read FASTA/FASTQ> [options]" << std::endl;
    std ::cout << " options:" << std::endl;
    std::cout << "    -m / --min-matches <min>   The minimum number of matches a read must have to pass [75]." << std::endl;
    std::cout << "    -t / --threads <threads>   The number of OpenMP threads to use [1]." << std::endl;
    std::cout << "    -s / --sketch-size <sketch size>   The size of the MinHash sketch to use per ref/read [1000]." << std::endl;
    std::cout << "    -k / --kmer-size <k>   The kmer size to use for hashing [16]" << std::endl;
    std::cout << "    -N / --ambiguous   Allow N's to be hashed along with canonical [ACTG] nucleotides [false]." << std::endl;
    std::cout << "    -z / --invert              Invert the query (i.e., return reads that *do not* match the references with min-matches) [false].." << std::endl;
    std::cout << "    -R / --ref-batch-size <batch>   The number of reference sequences to hash at one time (increase for speed at the cost of more memory usage) [10]." << std::endl;
    std::cout << "    -F / --fastq-batch-size <batch>   The number of read sequences to hash at one time (increase for speed at the cost of more memory usage) [1000]." << std::endl;
}

int main_filter(int argc, char** argv){
    std::vector<std::string> read_files;
    std::vector<std::string> ref_files;
    std::vector<std::string> pair_name_files;

    bool allow_ambiguous = false;
    int threads = 1;
    int kmer_length = 16;
    std::size_t min_length = 75;
    std::size_t sketch_size = 1000;
    std::size_t min_matches = 30;
    std::size_t ref_batch = 10;
    std::size_t read_batch = 1000;
    bool invert = false;
    // A hash_map to hold qname->read
    spp::sparse_hash_map<std::string, read_t> mate_cache;
    std::vector<read_t> ref_seqs;
    std::vector<hash_result> ref_hashes;

    optind = 2;
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"fastq", required_argument, 0, 'f'},
            {"min-matches", required_argument, 0, 'm'},
            {"reference", required_argument, 0, 'r'},
            {"ref-batch-size", required_argument, 0, 'R'},
            {"fastq-batch-size", required_argument, 0, 'F'},
            {"min-length", required_argument, 0, 'l'},
            {"pair-qnames", required_argument, 0, 'q'},
            {"sketch-size", required_argument, 0, 's'},
            {"ambiguous", no_argument, 0, 'N'},
            {"invert", no_argument, 0, 'z'},
            {"kmer-size", required_argument, 0, 'k'},
            {"version", no_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hzk:R:F:m:s:r:f:q:t:l:Nv", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'k':
                kmer_length = std::atoi(optarg);
                break;
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'R':
                ref_batch = std::atoi(optarg);
                break;
            case 'F':
                read_batch = std::atoi(optarg);
                break;
            case 'm':
                min_matches = std::atoi(optarg);
                break;
            case 'f':
                read_files.push_back(optarg);
                break;
            case 'q':
                pair_name_files.push_back(optarg);
                break;
            case 'z':
                invert = true;
                break;
            case 'N':
                allow_ambiguous = true;
                break;
            case 's':
                sketch_size = std::atoi(optarg);
                break;
            case 't':
                threads = std::atoi(optarg);
                break;
            case 'l':
                min_length = std::atoi(optarg);
                break;
            case '?':
            case 'h':
                usage(argv);
                exit(0);
            default:
                abort();
        }
    }

    omp_set_num_threads(threads);

    if (read_files.empty() & ref_files.empty()){
        std::cerr << "Please provide both a FASTA file of references and a FASTA/FASTQ file of reads." << std::endl;
        usage(argv);
    }

    /**
     * Read in an hash all the reference files.
     */
    //std::vector<read_t> curr_ref_fi_reads;
    //std::vector<hash_result> curr_ref_fi_hashes;
    for(std::size_t i = 0; i < ref_files.size(); ++i){
       
        klibpp::SeqStreamIn iss(ref_files[i].c_str());
        std::vector<klibpp::KSeq> records = iss.read(ref_batch);
        while (!records.empty()){
            #ifdef DEBUG_KRMR
            std::cerr << " records size: " << records.size() << std::endl;
            #endif
            std::size_t rsize = records.size();
            std::vector<read_t> curr_ref_fi_reads(rsize);
            std::vector<hash_result> curr_ref_fi_hashes(rsize);
            #pragma omp parallel for
            for (std::size_t i = 0; i < rsize; ++i){
                pliib::strcopy(records.at(i).name.c_str(), curr_ref_fi_reads.at(i).name);
                pliib::strcopy(records.at(i).seq.c_str(), curr_ref_fi_reads.at(i).seq);
                curr_ref_fi_reads[i].seqlen = records.at(i).seq.size();
                hash_sequence(curr_ref_fi_reads[i].seq,
                        curr_ref_fi_reads[i].seqlen,
                        kmer_length,
                        curr_ref_fi_hashes[i],
                        allow_ambiguous);
                curr_ref_fi_reads.at(i).purge_seq();            
                curr_ref_fi_hashes[i].sketch(sketch_size);
            }

            //#pragma omp critical
            {
                ref_seqs.insert(ref_seqs.end(), curr_ref_fi_reads.begin(), curr_ref_fi_reads.end());
                ref_hashes.insert(ref_hashes.end(), curr_ref_fi_hashes.begin(), curr_ref_fi_hashes.end());
            }
            records = iss.read(ref_batch);
        }
    }

    std::size_t n_ref_seqs = ref_seqs.size();
    std::cerr << "Processed " << ref_files.size() <<
        " files with " << n_ref_seqs <<
        " reference sequences." << std::endl;
    #ifdef DEBUG_KRMR
    #endif 

    std::size_t total_reads = 0;
    for(std::size_t i = 0; i < read_files.size(); ++i){
        klibpp::SeqStreamIn rss(read_files[i].c_str());
        std::vector<klibpp::KSeq> records = rss.read(read_batch);
        while (!records.empty()){
            std::size_t rsize = records.size();
            std::vector<read_t> curr_reads(rsize);
            std::vector<hash_result> curr_read_hashes(rsize);
            #pragma omp parallel for
            for (std::size_t i = 0; i < rsize; ++i){
                std::ostringstream st;
                char* seq = nullptr;
                std::size_t seq_len = records.at(i).seq.size();
                if (seq_len >= min_length){
                    pliib::strcopy(records.at(i).seq.c_str(), seq);
                    hash_sequence(seq,
                        seq_len,
                        kmer_length,
                        curr_read_hashes.at(i),
                        allow_ambiguous);
                    curr_read_hashes.at(i).sketch(sketch_size);
                    for (std::size_t j = 0; j < n_ref_seqs; ++j){
                        int inter = 0;
                        mkmh::hash_intersection_size(ref_hashes[j].hashes, ref_hashes[j].num_hashes,
                                curr_read_hashes.at(i).hashes, curr_read_hashes.at(i).num_hashes, inter);
                        //std::cerr << inter << std::endl;
                        if (invert ^ inter >= min_matches){
                            //std::cerr << "Record passed" << std::endl;
                            output_kseq(records.at(i), st);
                            st << std::endl;
                            //std::cout << records[i] << std::endl;;
                            std::cout << st.str();
                            st.str("");
                            break;
                        }
                    }
                    curr_read_hashes.at(i).clear();
                }
                else{
                    st << 
                        "Read filtered for length. QNAME: " << 
                        records.at(i).name << 
                        ", length: " << 
                        seq_len << std::endl;
                    std::cerr << st.str();
                    st.str("");
                }
                delete [] seq;
            }
            total_reads += rsize;
            std::cerr << "Processed " << total_reads << "  reads so far." << std::endl;
            records = rss.read(read_batch);
        }
    }

    std::cerr << "Processed " << total_reads << " total reads from " << read_files.size() << " files." << std::endl;


    return 0;

}
int main(int argc, char** argv){
    if (argc < 2){
        std::cerr << "No subcommand provided. Please provide a subcommand (try \" filter \")." << std::endl;
        usage(argv);
        return 1;
    }
     
    if (strcmp(argv[1], "filter") == 0){
        return main_filter(argc, argv);
    }
    else{
        std::cerr << "Invalid subcommand [" << argv[1] << "]. Please choose a valid subcommand (filter)" << std::endl;
        usage(argv);
    }
    return 0;
}
