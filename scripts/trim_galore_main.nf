

process trim_galore {

        cpus 16

        publishDir params.outdir2 , mode: 'copy'

        input:
		tuple val(sample), path(fastqs)

        output:
		path 'trim_galore_trimmed'

        conda "/etc/ace-data/home/hkawalya/miniconda3/envs/trim_galore"

        script:
        """
		trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired ${fastqs[0]} ${fastqs[1]} -o trim_galore_trimmed

        """
}
