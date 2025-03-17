//
// Alignment with Bwamem2
//

include { BWAMEM2_MEM         } from '../../modules/local/bwamem_align'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BWAMEM2 {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_index          // channel: /path/to/bwamem/index/
    save_unaligned    // val
    sort_bam          // val
    ch_fasta          // channel: /path/to/reference.fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    BWAMEM2_MEM ( ch_reads, ch_index, ch_fasta, save_unaligned, sort_bam )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( BWAMEM2_MEM.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig         = BWAMEM2_MEM.out.bam          // channel: [ val(meta), aligned ]
    log_out          = BWAMEM2_MEM.out.log          // channel: [ val(meta), log     ]
    fastq            = BWAMEM2_MEM.out.fastq        // channel: [ val(meta), fastq   ]

    bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi              = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                      // channel: [ versions.yml ]
}
