/*
========================================================================================
   Variant-Calling Nextflow Workflow
========================================================================================
   Description   : This pipeline performs fastqc and trimming with fastp
   Usage         : nextflow run fastp.nf
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
params.reads = "/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/*R{1,2}.fastq.gz"

println """\
         V A R I A N T-C A L L I N G - N F   P I P E L I N E
         ===================================================
         reads           : ${params.reads}
         outdir          : ${params.outdir}
         """
         .stripIndent()

/*
========================================================================================
   Create Channels
========================================================================================
*/

reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

    //FASTQC( reads_ch )
    Fastp ( reads_ch )
    // Enter the rest of the processes for variant calling based on the script

}

/*
========================================================================================
   Processes
========================================================================================
*/

/*
 * Perform a fastqc to guide your quality control process.

process FASTQC {
    tag{"FASTQC ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')

    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" ), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """
}

*/

//Trimming with fastp

//Defining the help function
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        "fastp -i seq.fastq -o seq_trimmed.fastq" or "fastp -i seq.R1.fastq -o seq_trimmed_R1.fastq.gz -I seq.R2.fastq -O seq_trimmed_R2.fastq.gz"

        Options:
        -i, --in1        read1 input file name (string [=])
        -o, --out1       read1 output file name (string [=])
        -I, --in2        read2 input file name (string [=])
        -O, --out2       read2 output file name (string [=])
  
       Optional arguments:
        -a, --adapter_sequence         The adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
            --adapter_sequence_r2      The adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])
            --adapter_fasta            Specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
            --detect_adapter_for_pe    By default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
        -V, --verbose                  Output verbose log information (i.e. when every 1M reads are processed).
        -l, --length_required          Reads shorter than length_required will be discarded, default is 15. (int [=15])
            --length_limit             Reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
        -f, --trim_front1              Trimming how many bases in front for read1, default is 0 (int [=0])
        -t, --trim_tail1               Trimming how many bases in tail for read1, default is 0 (int [=0])
        --help                         This usage statement.
        """
}

process Fastp {

        publishDir "${params.outdir}/Trimmed_reads_Fastp", mode: 'copy'

        input:
        tuple val( sample_id ), path( reads )

        output:
        path "*"

        conda 'fastp'
        script:
        """
         fastp \
                -i ${reads[0]} \
                -I ${reads[1]} \
                -o ${sample_id}_R1_trimmed.fastq.gz \
                -O ${sample_id}_R2_trimmed.fastq.gz
        """

}


workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
