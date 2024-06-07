nextflow.enable.dsl=2

process  fastqc {
	publishDir params.outdir , mode: 'copy'
	conda "$projectDir/conda_envs/env-preprocessing"
	input:
		path fastqs
	output:
		path '*', emit: qc_report
	script:
	"""
		fastqc $fastqs
	"""
}


process multiqc {
   publishDir params.outdir , mode: 'copy'
   conda "multiqc"
   input:
       path fastqcs
    output:
        path 'multiqc_report.html'
    script:
    """
    multiqc .
    """
}
