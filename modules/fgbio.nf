
process SET_MATE_INFO {
    tag "${subject}_${sample_id}"
    container ''
    input:
        tuple val(subject), val(sample_id),  path(bam), path(bai)
        path (reference)
    output:
        tuple val(subject), val(sample_id),
        path("${sample_id}.aligned.matefixed.bam")
    
    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir="/fs04/scratch2/vh83/jason/tmp" -jar "/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar" SortBam \
            -i $bam -o ${sample_id}.aligned.sorted.bam -s Queryname
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir="/fs04/scratch2/vh83/jason/tmp" -jar "/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar" SetMateInformation \
            -i ${sample_id}.aligned.sorted.bam -r $reference -o ${sample_id}.aligned.matefixed.bam
    """
}


process GROUP_READS {
    
    tag "${subject}_${sample_id}"
    container ''

    publishDir path: 'fix_this', mode: 'copy', pattern: "*.tsv"
    input:
        tuple val(subject), val(sample_id), path(bam)
    output:
        tuple val(subject), val(sample_id), path("${sample_id}.piped.grouped.histogram.tsv"), path("${sample_id}.piped.grouped.bam") 
    
    

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir="/fs04/scratch2/vh83/jason/tmp" -jar "/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar" GroupReadsByUmi \
        -i ${bam} -f "${sample_id}.piped.grouped.histogram.tsv" -o "${sample_id}.piped.grouped.bam" -s Adjacency -e 1 
    """

}

process GENERATE_CONSENSUS {
    
    tag "${subject}_${sample_id}"
    container ''
    input:
        tuple val(subject), val(sample_id), file(hist), file(bam)
    output:
        tuple val(subject), val(sample_id), file("${sample_id}.consensus.unmapped.bam") 
    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir="/fs04/scratch2/vh83/jason/tmp" -jar "/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar" CallMolecularConsensusReads \
        --input $bam --output ${sample_id}.consensus.unmapped.bam \
        --error-rate-post-umi 30 --min-reads 3
    """
}

process FGBIO_STATS {
    
    tag "${subject}_${sample_id}"
    container ''
    publishDir path: './output/metrics/fgbio', mode: 'copy'

    input:
        tuple val(subject), val(sample_id), file(bam) 
    output:
        tuple val(subject), val(sample_id), file("*") 
    
    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir="/fs04/scratch2/vh83/jason/tmp" -jar "/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar" CollectDuplexSeqMetrics \
        -i $bam -o ${sample_id} 
        
    """
}

process MAP_CONSENSUS {
    
    tag "${subject}_${sample_id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
    'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    publishDir path: './output/bams', mode: 'copy'
    input:
        tuple val(subject), val(sample_id), file(bam)
    output:
        tuple val(subject), val(sample_id), file("${sample_id}.consensus.aligned.bam") 


    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar "/fs04/vh83/local_software/picard.jar" SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR "/fs04/scratch2/vh83/jason/tmp" | \
    bwa mem -M -t ${task.cpus} -p $reference /dev/stdin > ${sample_id}.temp.bam
    samtools sort -o ${sample_id}.consensus.aligned.bam ${sample_id}.temp.bam
    """
}

process INDEX {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
    'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    publishDir path: './output/bams', mode: 'copy'
    
    input:
        tuple val(subject), val(sample_id), file(bam) 
    output:
        tuple val(subject), val(sample_id), file(bam), file("${sample_id}.consensus.aligned.bam.bai") 


    script:
    """
    samtools index $bam ${sample_id}.consensus.aligned.bam.bai
    """

}