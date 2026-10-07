// test run (outside of pipeline folder!):
// nextflow run Computational_Workflows_project_RoehreSchuele -profile test,docker --outdir test_results

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED }  from '../modules/nf-core/fastqc/main' // we include it again with an alias for the second run after trimming
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// Trimming:
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main' 
include { FASTP                  } from '../modules/nf-core/fastp/main'
// mark duplicate:
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'   
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main' 
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main' 

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cw-rnaseq-roehre-schuele_pipeline'

// MODULES NECESSARY FOR HISAT2 PREP - START (To be added into subworkflow once it runs and there is time)
include { HISAT2_EXTRACTSPLICESITES         } from '../modules/nf-core/hisat2/extractsplicesites' // → extracts splicing sites from a gtf file, Input for: HISAT Build /Align 
include { HISAT2_BUILD                      } from '../modules/nf-core/hisat2/build' // → Builds HISAT2 Index(set of files from reference genome for moreefficient aligning ), Creates input for HISAT2 ALign
include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../modules/nf-core/gunzip'
include {UNTAR as UNTAR_HISAT2_INDEX        } from '../modules/nf-core/untar'
//include { FASTQ_HISAT2_ALIGN                } from '../modules/nf-core/fastq_align_hisat2'
include { FASTQ_HISAT2_ALIGN                } from '../subworkflows/nf-core/fastq_align_hisat2'
// MODULES NECESSARY FOR HISAT2 PREP - END

//Salmon Quantification:
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../modules/nf-core/rsem/preparereference/main'
include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT } from '../modules/nf-core/salmon/quant/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CW_RNASEQ_ROEHRE_SCHUELE {

    take: // inputs for workflow
    ch_samplesheet // channel: samplesheet read in from --input
    //fasta // reference genome
    multiqc_config //could change report apperance in this file
    multiqc_logo // display multiqc logo in report
    multiqc_methods_description // optional text/template used for the methods section in the MultiQC report
    outdir
        // HISAT PREP VAR - START
    fasta
    gtf 
    //hisat2_index
    //splicesites
    //prepare_tool_indices // idk yet what this is for, probably for saving HISAT stuff
    // HISAT PREP VAR- END

    main:

    def ch_versions = channel.empty() // tracks the versions of tools
    def ch_multiqc_files = channel.empty() 

    // for new test with multiple fastq:
    CAT_FASTQ(ch_samplesheet)
    ch_merged_reads = CAT_FASTQ.out.reads

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_merged_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })


    //
    // MODULE: Trim reads with fastp (in larger pipelines we could wrap this in a subworkflow)
    //

    // first we need an input channel for fastp with all the required inputs
    ch_fastp_input = ch_merged_reads.map { meta, reads ->
        [meta, reads, []] 
        // fastp need meta data, reads and possible adapters as input
        // since we dont have adapter file we pass an empty list []
    }

    // then we can actually run fastp
    FASTP(
        ch_fastp_input,
        false,  // From documentation: Specify true to not write any reads that pass trimming thresholds. This can be used to use fastp for the output report only.
        false,  // Specify true to save files that failed to pass trimming thresholds ending in *.fail.fastq.gz
        false   // Specify true to save all merged reads to a file ending in *.merged.fastq.gz
    )
    // we can check which output fastp produces to know how to handle them:
    // sed -n '/output:/,/when:/p' modules/nf-core/fastp/main.nf

    // we also need to include version information and mix with multiQC!
    //ch_versions = ch_versions.mix(FASTP.out.versions_fastp) 

    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })
    //ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.map{ _meta, file -> file })

    // fastp output we want to pass the aligner:
    ch_trimmed_reads = FASTP.out.reads

    // rnaseq runs FASTQC here again and so do we:
    FASTQC_TRIMMED(ch_trimmed_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.map{ _meta, file -> file })
        
    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'cw-rnaseq-roehre-schuele_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    //
    // MODULE: HISAT2 BUILD AND SPLICESITE DETECTION
    //
    
    // Unzip fasta and gtf file using the processes created with GUNZIP Module
    // ONly allows GTF usage for now
    if (fasta.endsWith('.gz')) {
      ch_fasta = GUNZIP_FASTA ( [ [:], fasta]).gunzip.map {it[1]}
      //ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions) 
      } else {
      ch_fasta = Channel.value(file(fasta))}
    

    if (gtf.endsWith('.gz')) {
      ch_gtf = GUNZIP_GTF ( [ [:] , gtf ] ).gunzip.map { it[1] }
      //ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
      ch_gtf = Channel.value(file(gtf))
    }
  
    
    // Uncompress / Generate HISAT2 index from scratch
    

    ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
    //ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

    // Build index:
    ch_hisat2_index = HISAT2_BUILD ( ch_fasta
    .concat(ch_gtf)
    .concat(ch_splicesites) // put files in one stream
    .collect() // stream to list
    .map { fasta_file, gtf_file, splicesites_file ->
        [[:], fasta_file, gtf_file, splicesites_file]
    }, '7.GB').index.map { it[1] }
    //ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)


    ch_fasta_for_faidx = ch_fasta.map { fasta_file ->
        [
            [id: 'reference'],
            fasta_file,
            []
        ]
    }
    SAMTOOLS_FAIDX(
    ch_fasta_for_faidx,
    false
    )

    // combine fasta and fai as ecpected by hisat:
    ch_fasta_fai = ch_fasta
    .combine(
        SAMTOOLS_FAIDX.out.fai.map { meta, fai ->
            fai
        }
    )
    .map { fasta_file, fai_file ->
        [
            [id: 'reference'],
            fasta_file,
            fai_file
        ]
    }

// Add metadata back to HISAT2 index
ch_hisat2_index_for_align = ch_hisat2_index.map { index ->
    [[id: 'reference'], index]
}

// Add metadata back to splice sites
ch_splicesites_for_align = ch_splicesites.map { splicesites ->
    [[id: 'reference'], splicesites]
}

// Run alignment
FASTQ_HISAT2_ALIGN(
    ch_trimmed_reads,
    ch_hisat2_index_for_align,
    ch_splicesites_for_align,
    ch_fasta_fai,
    false
)

    //
    //MODULE: Quantify with Salmon
    //
    // salmon expects an transcript fasta file, in rnaseq this is build with RSM (we included RSM here as MAKE_TRANSCRIPTS_FASTA)
    ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta


    ch_salmon_index_input = ch_transcript_fasta
    .combine(ch_fasta)
    .map { transcript_fasta, genome_fasta ->
        [
            [id: 'reference'],
            transcript_fasta,
            genome_fasta
        ]
    }
    ch_salmon_index = SALMON_INDEX(ch_salmon_index_input).index
    
    //ch_salmon_index_for_quant = ch_salmon_index.map { index -> [[id: 'reference'], index]}

    // check which input is required:
    // sed -n '/input:/,/output:/p' \modules/nf-core/salmon/quant/main.nf input:
    ch_salmon_reference = ch_salmon_index
    .combine(ch_gtf)
    .combine(ch_transcript_fasta)
    .map { meta, index, gtf_file, transcript_fasta ->
        [
            meta,
            index,
            gtf_file,
            transcript_fasta
        ]
    }
    .first() // need this so not only the first file gets processed, need to figure out why this works

    SALMON_QUANT(
    ch_trimmed_reads, // tuple val(meta), path(reads)
    ch_salmon_reference //tuple val(meta2), path(index), path(gtf), path(transcript_fasta)
    )

    ch_multiqc_files = ch_multiqc_files.mix(
    SALMON_QUANT.out.json_info.map { meta, file -> file }
    )

    ch_multiqc_files = ch_multiqc_files.mix(
    SALMON_QUANT.out.lib_format_counts.map { meta, file -> file }
    )


    //
    //MODULE: Mark duplicate reads with picard_markduplicates
    //
    //could later add UMI step before this, enable with flag/parameter option
    // for picard we need bam (get this from Hisat, maybe need to sort it), fasta (reference, from parameter)
    //ch_sorted_bam = FASTQ_ALIGN_HISAT2.out.bam
    // and we need fai (createthis from ref fasta parameter with samtools) 
    // -> also get fai from Hisat subworkflow!
    // Hisat output: bam, log, fastq.gz,....
    // fasta = refernec genome
/*
    // get reference genome and prepare it for nfcore picard:
    ch_reference = Channel.value([
    [id: 'genome'],
    file(fasta)
    ])

    // use reference genome index for picard (fai file):
    SAMTOOLS_FAIDX(
    ch_reference,
    false
    )

    // sort bam for picard:
    SAMTOOLS_SORT(ch_bam) // here add Hisat output name!

    // run Picard markduplicates:
    PICARD_MARKDUPLICATES(
    ch_bam,
    ch_reference,
    SAMTOOLS_FAIDX.out.fai
    )
/*
    ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
    ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
    ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
    ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
    ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
    ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
    if (params.bam_csi_index) {
        ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
    }
*/

    //
    // MODULE: MultiQC (summary step at the end!)
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'cw-rnaseq-roehre-schuele'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    //versions       = ch_versions  // channel: [ path(versions.yml) ]
    fasta = ch_fasta // channel: path(genome.fasta)
    gtf = ch_gtf // channel: path(genome.gtf)
    splicesites = ch_splicesites //channel: path(genome.splicesites)
    hisat2_index = ch_hisat2_index // channel: path(hisat2/index/)
    salmon_results = SALMON_QUANT.out.results

    
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
