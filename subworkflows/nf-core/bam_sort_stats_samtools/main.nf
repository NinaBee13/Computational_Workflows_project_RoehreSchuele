//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]

    main:
    // we need to make the reference reusable for every BAM
    // we take first item from queue channel (item = whole tuple [val(meta), path(fasta), path(fai) ])
    // and .first() transforms this into a reusbale value channel
    ch_fasta_fai_value = ch_fasta_fai.first()

    SAMTOOLS_SORT(
        ch_bam, 
        //ch_fasta_fai,
        ch_fasta_fai_value,
         ''
    )

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.index, by: [0])
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS(
        ch_bam_bai, 
        //ch_fasta_fai
        ch_fasta_fai_value
    )

    emit:
    bam      = SAMTOOLS_SORT.out.bam // channel: [ val(meta), [ bam ] ]
    index    = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
