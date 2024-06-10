nextflow.enable.dsl=2

//Defining the help function

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --PE or --SE "Path to the reads"  < Other options >

        Mandatory arguments:
         --PE	                        paired end reads (fastq files)
         --SE				Single end reads (fastq file)
         --ref                          Reference genome (full path required)
         --adapter                      Adapter sequence to use for trimming (required for trimmomatic tool)

       Optional arguments:
        --outdir                       Output directory to place final BLAST output
        --phred			       Phred score (default = phred 33) for trimmomatic
        --phred2		       Phred score (default = Phred 33) for trim_galore
        --options                      Adapter sequence for trim_galore [default illumina]
        --threads                      Number of CPUs to use during blast job [16]
        --tool                         Which tool to use to trim reads [trimmomatic, trimgalore, fastp]
        --LEADING		       Cut bases off the start of a read, if below a threshold quality
        --TRAILING		       Cut bases off the end of a read, if below a threshold quality
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



      
// Initial fastqc 

process Fastqc_PE {
    //cpus 4
    //memory '6 GB'
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
process Fastqc_SE {
    //cpus 4
    //memory '6 GB'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("*_fastqc.html")
    tuple val(sample), path("*_fastqc.zip")
    // path "multiqc_report.html"
    
    conda 'fastqc multiqc'
  

    script:
    """
    fastqc -o . ${reads} 
    # multiqc .
    """
}


// Trimming reads using trimmomatic

process Trimmomatic_PE {

    //cpus 4
    //memory '6 GB'
    publishDir "${params.outdir}/Trimmed_reads", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_1.trimmed.fastq.gz") 
    tuple val(sample), path("*_2.trimmed.fastq.gz") 
  
 
    conda 'trimmomatic'

    script:
    """
  
        trimmomatic  PE -threads $params.threads $params.phred ${reads[0]} ${reads[1]} \
        ${sample}_1.trimmed.fastq.gz ${sample}_1up.trimmed.fastq.gz ${sample}_2.trimmed.fastq.gz ${sample}_2up.trimmed.fastq.gz \
        ILLUMINACLIP:$params.adapter:2:30:10:2:True $params.LEADING $params.TRAILING $params.SLIDINGWINDOW $params.MINLEN
        
        rm -rf ${sample}_1up.trimmed.fastq.gz ${sample}_2up.trimmed.fastq.gz
           
    """
}

process Trimmomatic_SE {

    //cpus 4
    //memory '6 GB'
    publishDir "${params.outdir}/Trimmed_reads", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.trimmed.fastq.gz") 
  
 
    conda 'trimmomatic'

    script:
    """
  
        trimmomatic  SE -threads $params.threads $params.phred ${reads}  \
        ${sample}.trimmed.fastq.gz  \
        ILLUMINACLIP:$params.adapter:2:30:10:2:True \
        $params.LEADING $params.TRAILING $params.SLIDINGWINDOW $params.MINLEN
        
        rm -rf ${sample}_1up.trim.fastq.gz 
           
    """
}

// Trimming with trim-galore for both paired and single ended reads

process Trimgalore_PE {
        
        publishDir "${params.outdir}/Trimmed_reads_Trimgalore", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_1_val_1.fq.gz") 
	tuple val (sample),path("*_2_val_2.fq.gz") 
	
	conda 'trim-galore'
	script:
	"""
	
	trim_galore --paired  \
	 --three_prime_clip_R1 $params.INT \
	 --three_prime_clip_R2 $params.INT \
	 ${reads[0]} ${reads[1]} --fastqc
	
	"""

}

process Trimgalore_SE {
        
        publishDir "${params.outdir}/Trimmed_reads_Trimgalore", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	tuple val(sample), path ("*trimmed.fq.gz")
	
	conda 'trim-galore'
	script:
	"""
	
	trim_galore ${reads}  --basename ${sample}_trimmed.fq.gz --fastqc
	
	"""

}

//Trimming with fastp for both paired end and single reads

process Fastp_PE {
        
        publishDir "${params.outdir}/Trimmed_reads", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_1.trimmed.fastq.gz")
        tuple val(sample), path("*_2.trimmed.fastq.gz") 
	
	conda 'fastp'
	script:
	"""
	
	fastp -i ${reads[0]} \
	-o ${sample}_1.trimmed.fastq.gz  \
	-I ${reads[1]} \
	-O ${sample}_2.trimmed.fastq.gz \
	--detect_adapter_for_pe
	"""

}

process Fastp_SE {
        
        publishDir "${params.outdir}/Trimmed_reads", mode: 'copy'
        
	input:
        tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_1.trimmed.fastq.gz")
	
	conda 'fastp'
	script:
	"""
	
	fastp -i ${reads} -o ${sample}_1.trimmed.fastq.gz 
	"""

}

// Quality control after trimming for trimmomatic and fastp

process Fastqc2_PE {
    
    publishDir "${params.outdir}/fastqc_post", mode: 'copy'
    
    input:
    tuple val(sample), path(trimmed_R1)
    tuple val (sample),path(trimmed_R2) 
  
    
    output:
     
    tuple val(sample), path("*.trimmed_fastqc.html")
    tuple val(sample), path("*.trimmed_fastqc.zip")
    path "multiqc_report.html"  
    
    conda 'fastqc multiqc'


    script:
    """
    fastqc -o . ${trimmed_R1} ${trimmed_R2}
    multiqc .
    """
}

process Fastqc2_SE {
    
    publishDir "${params.outdir}/fastqc_post", mode: 'copy'
    
    input:
    tuple val(sample), path(trimmed_R1)
  
    
    output:
     
    tuple val(sample), path("*.trimmed_fastqc.html")
    tuple val(sample), path("*.trimmed_fastqc.zip")
    path "multiqc_report.html"  
    
    conda 'fastqc multiqc'


    script:
    """
    fastqc -o . ${trimmed_R1} 
    multiqc .
    """
}

// Checking if indexing the reference genome is required

//Indexing the reference genome

process index {

	publishDir "${params.index_dir}", mode: 'copy'
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

//Read mapping and conversion to bam

process alignment_bam_PE {
	publishDir "${params.outdir}/Bam", mode: 'copy'

	input:
	
	path index_dir
	path ref
	tuple val(sample), path(trimmed_R1)
	tuple val(sample),  path(trimmed_R2)	
	
	output:
	tuple val(sample), path ("*.bam")
	
	conda 'bwa'
	
	script:
	
	"""
	
	bwa mem ${index_dir}/${ref} ${trimmed_R1} ${trimmed_R2} | samtools view -S -b  -o ${sample}.bam -
	
		
	"""
}

//Read mapping and conversion to bam

process alignment_bam_SE {
	publishDir "${params.outdir}/Bam", mode: 'copy'

	input:
	
	path index_dir
	path ref
	tuple val(sample),path(trimmed_R1)	
	
	output:
	tuple val(sample), path ("*.bam")
	
	conda 'bwa'
	
	script:
	
	"""
	
	bwa mem ${index_dir}/${ref} ${trimmed_R1} | samtools view -S -b  -o ${sample}.bam -
	
		
	"""
}

// Converting the bam to sorted bam

process alignment_sorted_bam {
	publishDir "${params.outdir}/Sorted_bam", mode: 'copy'

	input:
	tuple val(sample), path (Bam)
	
	output:
	tuple val(sample), path ("*.sorted.bam")
	
	
	conda 'samtools'
	
	script:
	"""
	
	samtools sort -o ${sample}.sorted.bam ${Bam} 
	
	
	"""
}

// Checking the alignment stats of the sorted bam

process alignment_statistics {
	publishDir "${params.outdir}/Alignment_Stats", mode: 'copy'

	input:
	tuple val(sample), path (Sorted_bam)

	
	output:
	tuple val(sample), path ("*_stats.txt")
	
	
	conda 'samtools'
	
	script:
	"""
	
	
	
	samtools flagstat ${Sorted_bam} > ${sample}_stats.txt
	
		
	"""
}


workflow {

    // Define input channels 
    
    
    
    read_pairs_ch = Channel.fromFilePairs(params.Paired_end_reads, checkIfExists: true)
                            .ifEmpty { error "Cannot find any FASTQ files in: ${params.Paired_end_reads}" }
    
    // Define the channel for single end reads, Maps each file in the channel to a tuple containing the base name of the file (sample name) and the file itself 
    // Removes suffix like "_1.fastq.gz"
    
    			   
    read_SE_ch = Channel.fromPath(params.Single_end_reads, checkIfExists: true)
                   .ifEmpty { error "Cannot find any FASTQ files in: ${params.Single_end_reads}" }
                   .map { file -> 
                       def fileName = file.getName()
                       def sample = fileName.replaceAll(/_1\.fastq\.gz/, "")
                       tuple(sample, file)
                   }			   

    ref_ch = Channel.fromPath(params.ref)
    
    //gbk_ref_ch = Channel.fromPath(params.gbk_ref)
    
    index_ch=Channel.fromPath(params.index_dir)
  
    
    if (params.PE) {
        
        Fastqc_PE(read_pairs_ch)
        
    } else if (params.SE) {
    
	Fastqc_SE(read_SE_ch)
    }
    
    
    // Choose  the type of reads as either paired end or single end reads and also which tool to use for trimming 
    
    if (params.PE && params.tool == 'trimmomatic') {
        trim_ch = Trimmomatic_PE(read_pairs_ch)
        
    } else if (params.SE && params.tool == 'trimmomatic'){
    	trim_ch = Trimmomatic_SE(read_SE_ch)
             
    } else if (params.PE && params.tool == 'trim_galore') {
    	trim_ch = Trimgalore_PE(read_pairs_ch)
    
    } else if (params.SE && params.tool == 'trim_galore') {
    	trim_ch = Trimgalore_SE(read_SE_ch)
    
    } else if (params.PE && params.tool == 'fastp') {
    	trim_ch = Fastp_PE(read_pairs_ch)
    
    } else if (params.SE && params.tool == 'fastp') {
        trim_ch = Fastp_SE(read_SE_ch)
        
    } 

    // Index reference genome
    
    index(ref_ch)
    
    // Optional step Perform FastQC for Trimmomatic and Fastp
      
    if (params.PE && params.tool == 'trimmomatic') {
        Fastqc2_PE(Trimmomatic_PE.out)
        
    } else if (params.PE && params.tool == 'trimmomatic') {
	Fastqc2_SE(Trimmomatic_SE.out)
		
    } else if (params.PE && params.tool == 'fastp') {
	Fastqc2_PE(Fastp_PE.out)	
	
    }else if (params.SE && params.tool == 'fastp') {
	Fastqc2_SE(Fastp_SE.out)	
    }
    
    
    //alignment_bam(index_ch, ref_ch, Fastp.out)
    
    if (params.PE) {
        
        bam_ch = alignment_bam_PE(index_ch, ref_ch, trim_ch)
        
    } else if (params.SE) {
    
	bam_ch = alignment_bam_SE(index_ch, ref_ch, trim_ch)
    }
    
    
    sorted_ch = alignment_sorted_bam(bam_ch)
    
    stats_ch = alignment_statistics(sorted_ch)
    
    
}
