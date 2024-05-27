nextflow.enable.dsl=2

//Defining the help function

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "Path to the reads"  --ref "reference genome" 

        Mandatory arguments:
         --reads                        paired end reads (fastq files)
         --ref                          Reference genome (full path required)
         --adapter                      Adapter sequence to use for trimming (required for trimmomatic tool)

       Optional arguments:
        --outdir                       Output directory to place final output
        --phred			       Phred score (default = phred 33) for trimmomatic
        --phred2		       Phred score (default = Phred 33) for trim_galore
        --options                      Adapter sequence for trim_galore [default illumina]
        --threads                      Number of CPUs to use during blast job [16]
        --tool                         Which tool to use to trim reads [trimmomatic, trimgalore, fastp]
        --LEADING		       Cut bases off the start of a read, if below a threshold quality
        --TRAILING		       Cut bases off the end of a read, if below a threshold quality
        --HEADCROP		       Cut the specified number of bases from the start of the read
        --SLIDINGWINDOW		       Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls belo            
                                       below a threshold [deafult SLIDINGWINDOW:4:15]
        --MINLEN		       Drop the read if it is below a specified length 
        --help                         This usage statement.
        """
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}



      
///fastqc 

process Fastqc {
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("*_fastqc.html")
    tuple val(sample), path("*_fastqc.zip")
    path "multiqc_report.html"
    
    conda 'fastqc multiqc'
  

    script:
    """
    fastqc -o . ${reads[0]} ${reads[1]}
    multiqc .
    """
}


// Trimming reads using trimmomatic

process Trimmomatic {

    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/Trimmed_reads_Trimmomatic", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_1.trim.fastq.gz")
    tuple val(sample), path("*_2.trim.fastq.gz")
  
 
    conda 'trimmomatic'

    script:
    """
  
        trimmomatic  PE -threads $params.threads $params.phred ${reads[0]} ${reads[1]} \
        ${sample}_1.trim.fastq.gz ${sample}_1up.trim.fastq.gz ${sample}_2.trim.fastq.gz ${sample}_2up.trim.fastq.gz \
        ILLUMINACLIP:$params.adapter:2:30:10:2:True $params.LEADING $params.TRAILING $params.SLIDINGWINDOW $params.MINLEN
        
        rm -rf ${sample}_1up.trim.fastq.gz ${sample}_2up.trim.fastq.gz
           
    """
}



// Quality control after trimming for trimmomatic 

process Fastqc2 {
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/fastqc_post", mode: 'copy'
    
    input:
    tuple val(sample), path("*_1.trim.fastq.gz")
    tuple val(sample), path("*_2.trim.fastq.gz")
    
    output:
    tuple val(sample), path("*_fastqc.html")
    tuple val(sample), path("*_fastqc.zip")
    path "multiqc_report.html"
    
    
    // conda '/home/sandra/miniconda3/envs/TBprofiler2'

    script:
    """
    fastqc -o . ${reads[0]} ${reads[1]}
    multiqc .
    """
}
// Trimming with trim-galore

process Trimgalore {
        
        publishDir "${params.outdir}/Trimmed_reads_Trimgalore", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	path "*"
	
	conda 'trim-galore'
	script:
	"""
	
	trim_galore --paired   --three_prime_clip_R1 $params.INT --three_prime_clip_R2 $params.INT ${reads[0]} ${reads[1]} --fastqc
	
	"""

}

//Trimming with fastp
process Fastp {
        
        publishDir "${params.outdir}/Trimmed_reads_Fastp", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	path "*"
	
	conda 'fastp'
	script:
	"""
	
	fastp -i ${reads[0]} -o ${sample}_1.trim.fastq.gz  -I ${reads[1]} -O ${sample}_2.trim.fastq.gz
	
	"""

}


//Indexing the reference genome

process index {


	input:
	path ref_genome
	
	output:
	path "*"
	
	conda 'bwa'
	
	script:
	"""
	bwa index ${ref_genome} 
		
	"""
}


workflow {
    // Define input channels
    Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any FASTQ files in: ${params.reads}" }
        .set { read_pairs_ch }
      
    Channel.fromPath(params.ref).set {ref_ch}
    Fastqc(read_pairs_ch)
    
    if (params.tool == 'trimmomatic') {
 		Trimmomatic(read_pairs_ch)    
    } 
    else if (params.tool == 'trim_galore') {
		Trimgalore(read_pairs_ch)
    
    }
    else{
		Fastp(read_pairs_ch)
    }
    index(ref_ch)
   
   
}
