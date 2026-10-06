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
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cw-rnaseq-roehre-schuele_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CW_RNASEQ_ROEHRE_SCHUELE {

    take: // inputs for workflow
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config //could change report apperance in this file
    multiqc_logo // display multiqc logo in report
    multiqc_methods_description // optional text/template used for the methods section in the MultiQC report
    outdir

    main:

    def ch_versions = channel.empty() // tracks the versions of tools
    def ch_multiqc_files = channel.empty() 
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })


    //
    // MODULE: Trim reads with fastp (in larger pipelines we could wrap this in a subworkflow)
    //

    // first we need an input channel for fastp with all the required inputs
    ch_fastp_input = ch_samplesheet.map { meta, reads ->
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
    //ch_versions = ch_versions.mix(FASTP.out.versions_fastp) // I think we dont need to collect ther version information ourselfs bc of the topic channel!

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
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
