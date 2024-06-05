nextflow.enable.dsl=2

params.fastqcs = '/etc/ace-data/CancerGenomicsWG/users/mnamuswe/results/fastqc_reports'
params.outdir  = '/etc/ace-data/CancerGenomicsWG/users/mnamuswe/results/fastqc_reports'

include { multiqc } from './main.nf'

workflow {
	fastqc_s = channel.fromPath( params.fastqcs )

	multiqc( fastqc_s )
}
