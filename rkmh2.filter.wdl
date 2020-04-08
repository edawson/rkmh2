task rkmhFilterReads{
   input{
    Array[File] refs
    Array[File] reads
    Int threads
    Int diskGB
    }
    command {

    }
    runtime {

    }
    output{

    }
}

workflow rkmh2Filter{
    Array[File] refs
    Array[File] reads
    Int? threads = 8
    Int diskGB = ceil(size(refs, "gb") + size(reads, "gb") * 2)

    call rkmhFilterReads as filterReads{
        input:
            refs=refs,
            reads=reads,
            threads=threads,
            diskGb=diskGB
    }
}
