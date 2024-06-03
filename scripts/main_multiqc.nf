nextflow.enable.dsl=2


process multiqc {

   publishDir params.outdir , mode: 'copy'

   input:
       path fastqcs 
 
    output:
        path 'multiqc_report.html'
  
 
    script:
    """
    multiqc .
    """
}
