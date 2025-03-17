include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_MULTIMAP } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MULTIMAP } from '../../modules/nf-core/samtools/index/main'
include { FILTER_LOW_COV } from '../../modules/local/filter_low_cov'
include { SUBREAD_FEATURECOUNTS } from '../../modules/nf-core/subread/featurecounts/main'
include { MERGE_COUNTS } from '../../modules/local/merge_counts'

workflow VIRAL_READ_COUNT {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai ] ]
    fasta   // path   : genome.fasta
    gff     // path   : genome.gff

    main:

    ch_versions = Channel.empty()

    //
    // Filter multi-mapped reads using Samtools view
    //
    SAMTOOLS_FILTER_MULTIMAP (
        bam_bai,
        fasta,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FILTER_MULTIMAP.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX_MULTIMAP (
        SAMTOOLS_FILTER_MULTIMAP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MULTIMAP.out.versions.first())

    SAMTOOLS_FILTER_MULTIMAP.out.bam
    | join(SAMTOOLS_INDEX_MULTIMAP.out.bai, by: [0])
    | set { ch_bam_bai }

    FILTER_LOW_COV (
        ch_bam_bai,
        fasta
    )
    ch_versions = ch_versions.mix(FILTER_LOW_COV.out.versions)

    FILTER_LOW_COV.out.high_coverage_bam
    | combine (gff)
    // | map { meta, bam, gff -> tuple(meta, bam, gff) }
    | set { ch_bam_annotation }

    SUBREAD_FEATURECOUNTS (
        ch_bam_annotation
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    SUBREAD_FEATURECOUNTS.out.counts
    | map { meta, count -> count }
    | toList()
    | set { ch_counts }

    MERGE_COUNTS (
        ch_counts
    )

    emit:
    merged_counts = MERGE_COUNTS.out.merged_counts       // path: merged_counts.csv
    bam    = FILTER_LOW_COV.out.high_coverage_bam  // channel: [ meta, *.bam ]
    versions = ch_versions                          // channel: [ versions.yml ]
}
