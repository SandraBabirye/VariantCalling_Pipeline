nextflow.enable.dsl=2

params.R1 = '/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/*R1.fastq.gz'
params.R2 = '/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/*R2.fastq.gz'
params.ref = '/etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/Homo_sapiens_assembly38.fasta'
params.index = '/etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/Homo_sapiens_assembly38.fasta*'
params.outdir = '/etc/ace-data/CancerGenomicsWG/users/gnakabiri/results'

include { bwa_index } from './RGmain.nf'
include { bwa_mem } from './RGmain.nf'

workflow {
    ref = channel.fromPath(params.ref)
    indices = channel.fromPath(params.index).collect()
    fastq_R1 = Channel.fromPath(params.R1)
    fastq_R2 = Channel.fromPath(params.R2)
    sample_pairs = fastq_R1.map { path ->
        def sample_id = path.getName().tokenize('_')[0]  // Assuming sample ID is part of the filename before the first underscore
        def fastq_R2 = path.toString().replace('R1.fastq.gz', 'R2.fastq.gz')
        return tuple(sample_id, path, file(fastq_R2))
     }

     bwa_mem(sample_pairs, ref, indices)
}
    
