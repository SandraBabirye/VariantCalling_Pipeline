

include { trim_galore } from './trim_galore_main.nf'

workflow {
    reads = Channel.fromFilePairs(params.fastqs)

    trim_galore(reads)
}
