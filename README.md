rkmh2: Kmer-based read filtering using MinHash
----------------------------------------------

![C/C++ build/test for rkmh2](https://github.com/edawson/rkmh2/workflows/C/C++%20build/test%20for%20kramer/badge.svg)
![docker: hpobiolab/rkmh2](https://img.shields.io/docker/image-size/hpobiolab/rkmh2)

## Introduction
**rkmh2** is a reimplementation of (some) of the algorithms for MinHash-based read filtering from the
original [rkmh package](https://github.com/edawson/rkmh) and [accompanying paper](https://doi.org/10.1186/s12859-019-2918-y). It focuses on improved stability, lower memory usage, usability and speed.

### Improvements over the original rkmh
- [x] batched ref / read hashing (greatly reduces memory usage)
- [x] better code flow / fewer branches (speed - more gains to be made, too)
- [x] safer code (I learned more C++)
- [x] simpler code - we're down to essentially four header files and one cpp.
- [x] big improvements in memory usage (WGS human, 30k hashes per sequence on standard laptops)
- [x] Native fastq output for filter without ever saving the sequences.

### Notably missing
- [ ] no `rkmh classify` or `stream` replacements
- [ ] no `-I/-M` depth-based filters (which were
    the major improvements of rkmh over other methods to begin with)
- [ ] no streaming input of reads yet
- [ ] code not yet profiled
- [ ] no HPV16-specific features

### Features to come
- [ ] paired-end read hashing / filtering
- [ ] replace the original `rkmh classify` command
- [ ] further speed improvements, new hashing algos, and more parallelization
- [ ] automated tests (in the works) and a cleaner repository.


## filter : select/filter reads matching a set of sequences
Currently only one command is supported: `rkmh2 filter`

### filter basics
*rkmh filter* performs "select" or "filter" operations on a readset, depending on options. It requires
at a minimum a set of FASTA/FASTQ references (which may be reads) and a set of sequences to filter.

The default behavior is to return reads which passed filtering criteria. The default filter is
a minimum of 75 matches using a 1000 hash MinHash sketch. The sketch size can be adjusted using the `-s <SKETCH>` 
parameter and the minimum required number of matches can be adjusted with `-m <MINMATCH>`. 
Here's an example that uses a sketch size of 4000 (~40% of the Zika ref genome length) and 
a minimum of 15 matches:
```
./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -m 10 -s 4000 > matches.fq
```

The output of `rkmh2 filter` is again a FASTQf file, so if you wanted to filter the prefiltered matches.fq file to use a more stringent minimum of
15 matches you could:
```
./rkmh2 filter -r data/zika.refs.fa -f matches.fq -m 15 -s 4000 > strict_matches.fq
```

### Length filtering
You can also remove reads shorter than a certain length using `-l <MINLENGTH>. If you wanted to remove
reads shorter than 50bp, for example:
```
 time ./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -l 50 -m 0 > matches.fq
```
Just make sure to set the `-m <MINMATCH>` parameter to zero if you want to filter exclusively
on the length of a read.

This can be combined with sketch filtering or used independently. Length filtering cannot
yet utilize the `-z` flag to return reads *shorter* than a given length, but that *could* be a future feature.

### Performance tuning options : sketch size, threads, batch size
Increasing the sketch size increases sensitivity for long sequences, at the cost of doing lots more work.
At the extreme, a very long sketch size will devolve into just exact kmer matching (where every kmer in a read is compared against
every kmer in a reference sequence).

You can use large sketch sizes for large references - tens of thousands of hashes may be in a sketch. However,
run time will scale linearly with sketch size.
```
 time ./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -s 30000 > matches.fq
```

rkmh2 is multithreaded with OpenMP, with plans to expand into other parallelization technologies. Use the `-t <THREADS>`
option to use multiple threads:
```
 time ./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -m 15 -s 4000 -t 4 > matches.fq
```

The number of sequences processed per batch is also adjustable. Low batch sizes save significantly on memory use at the cost of some overhead
caused by having to read off disk and start the OpenMP loop. `-R <REFBATCH>` adjusts the number of references hashed at once;
`-F <READBATCH>` adjusts the number of reads hashed at once.
```
 time ./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -m 15 -s 4000 -R 100 -F 100000 > matches.fq
```

Defaults are listed in the usage, which is accessible by:
```
./rkmh2 filter -h
```

### Inverting searches
The `-z` flag will return sequences which do not meet the minimum matches requirements. Think of this flag like grep's `-v`:
```
 time ./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -s 4000 -m 15 -t 4 -z > non-matches.fq
```

## Future command development
My hope is to revive the original rkmh classify command as well as implement the kmer-counting filters applied in the paper.


