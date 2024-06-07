nextflow.enable.dsl=2

params.fastqs = '/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/DR-077-20_R*fastq.gz'
params.outdir  = '/etc/ace-data/CancerGenomicsWG/VariantCalling/pipeline/Results'

include { fastqc } from './main.nf'
include { multiqc } from './main.nf'

workflow {
	fastq = channel.fromPath( params.fastqs )
	fastqc( fastq )
	multiqc(fastqc.out.qc_report)
}
