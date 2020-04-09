version 1.0

task rkmhFilterReads{
   input{
    File refs
    File reads
    Int? sketchSize = 1000
    Int? minMatches = 10
    Boolean? invert = false
    Int? minLength = 50
    Int? refBatch
    Int? readBatch
    Int? threads = 1
    Int diskGB
    String outbase = basename(basename(basename(basename(reads, ".gz"), ".fq"), "fastq"), ".fa")
    }
    command {
      rkmh2 ~{true="-z" false="" invert} \
        -f ~{reads} \
        -r ~{refs} \
        ~{"-R " + refBatch} ~{"-F " + readBatch} -t ~{threads} \
        ~{"-l " + minLength} ~{"-s " + sketchSize} ~{"-m " + minMatches} > ~{outbase}.rkmhFiltered.fastq
    }
    runtime {
      docker : "hpobiolab/rkmh2"
      cpu : threads
      memory : "16GB"
      disks : "local-disk " + diskGB + " HDD"
      preemptible : 3
    }
    output{
     File filteredReads = "~{outbase}.rkmhFiltered.fastq" 
    }
}

workflow rkmh2Filter{
  input{
    File refs
    File reads
    Int? refBatch
    Int? readBatch
    Int? sketchSize
    Int? minMatches
    Boolean? invert
    Int? threads = 8
    Int diskGB = ceil(size(refs, "gb") + size(reads, "gb") * 2)
  }

    call rkmhFilterReads as filterReads{
        input:
            refs=refs,
            reads=reads,
            refBatch=refBatch,
            readBatch=readBatch,
            sketchSize=sketchSize,
            minMatches=minMatches,
            invert=invert,
            threads=threads,
            diskGB=diskGB
    }
}
