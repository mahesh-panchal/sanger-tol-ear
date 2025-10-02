/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subpipeline imports
// include { SANGER_TOL_BTK            } from '../modules/local/sanger-tol/blobtoolkit/main'
// include { SANGER_TOL_CPRETEXT       } from '../modules/local/sanger-tol/curationpretext/main'
include { BTK_INPUT } from '../modules/local/btk/input/main'
include { CPRETEXT_INPUT } from '../modules/local/cpretext/input/main'
include { NEXTFLOW_RUN as SANGER_TOL_BTK } from '../modules/local/nextflow/run/main'
include { NEXTFLOW_RUN as SANGER_TOL_CPRETEXT } from '../modules/local/nextflow/run/main'

// Module imports
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
include { GENERATE_SAMPLESHEET } from '../modules/local/generate_samplesheet/main'
include { GFASTATS } from '../modules/nf-core/gfastats/main'
include { MERQURYFK_MERQURYFK } from '../modules/nf-core/merquryfk/merquryfk/main'

// Plugin imports
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ear_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EAR {
    take:
    ch_sample_id
    ch_reference_hap1
    ch_reference_hap2
    ch_reference_haplotigs
    ch_fastk_hist
    ch_fastk_ktab
    ch_longread_dir
    ch_cpretext_hic_dir
    ch_cpretext_telomotif
    ch_cpretext_aligner
    ch_btk_read_layout
    ch_btk_un_diamond_db
    ch_btk_nt_db
    ch_btk_ncbi_taxonomy_path
    ch_btk_taxid
    ch_busco_lineages
    ch_busco_config

    main:
    ch_versions = Channel.empty()

    //
    // NOTE: THIS STAYS HERE | MOVING IT INTO PIPELINE INIT BREAKS IT
    // LOGIC: SPLITS INPUT STEPS INTO A LIST THAT CONTROLS PROCESSES ON EXISTENCE
    //
    exclude_steps = params.steps ? params.steps.tokenize(",") : "NONE"
    full_list = ["btk", "cpretext", "merquryfk", "NONE"]

    if (!full_list.containsAll(exclude_steps)) {
        error("There is an extra argument given on Command Line: \nCheck contents of: ${exclude_steps}\nMaster list is: ${full_list}")
    }

    //
    // LOGIC: IF HAPLOTIGS IS EMPTY THEN PASS ON HALPLOTYPE ASSEMBLY
    //          IF HAPLOTIGS EXISTS THEN MERGE WITH HAPLOTYPE ASSEMBLY
    //
    ch_reference_haplotigs
        .ifEmpty('NO_HAPLOTIGS')  // Use a marker value for empty case
        .combine(ch_sample_id)
        .combine(ch_reference_hap2)
        .branch { haplotigs, sample_id, hap2 ->
            concat_needed: haplotigs != 'NO_HAPLOTIGS'
                return tuple([id: sample_id], [hap2, haplotigs])
            no_concat: true
                return hap2
        }
        .set { processing_branch }

    CAT_CAT(processing_branch.concat_needed)
    ch_versions = ch_versions.mix(CAT_CAT.out.versions)

    ch_haplotype_fasta = CAT_CAT.out.file_out
        .mix(processing_branch.no_concat)

    //
    // MODULE: ASSEMBLY STATISTICS FOR THE FASTA
    //
    GFASTATS(
        ch_reference_hap1,
        "fasta",
        [],
        [],
        [[], []],
        [[], []],
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions)


    //
    // LOGIC: STEP TO STOP MERQURY_FK RUNNING IF SPECIFIED BY USER
    //
    if (!exclude_steps.contains("merquryfk")) {
        //
        // LOGIC:  REFORMAT A BUNCH OF CHANNELS FOR MERQUERYFK
        //
        ch_reference_hap1
            .combine(ch_haplotype_fasta)
            .combine(ch_fastk_hist)
            .combine(ch_fastk_ktab)
            .map { meta1, primary, meta2, haplotigs, fastk_hist, fastk_ktab ->
                tuple(
                    meta1,
                    fastk_hist,
                    fastk_ktab,
                    primary,
                    haplotigs,
                )
            }
            .set { merquryfk_input }

        //
        // MODULE: MERQURYFK PLOTS OF GENOME
        //
        MERQURYFK_MERQURYFK(
            merquryfk_input,
            [],
            [],
        )
        ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions)
    }


    //
    // LOGIC: STEP TO STOP BTK RUNNING IF SPECIFIED BY USER
    //
    if (!exclude_steps.contains("btk")) {
        //
        // MODULE: GENERATE_SAMPLESHEET creates a csv for the blobtoolkit pipeline
        //
        GENERATE_SAMPLESHEET(
            ch_reference_hap1,
            ch_longread_dir,
            ch_btk_read_layout,
        )
        ch_versions = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

        //
        // MODULE: Run Sanger-ToL/BlobToolKit
        //
        SANGER_TOL_BTK(
            'sanger-tol/blobtoolkit',
            [
                "-profile ${workflow.profile}",
                "-r ${params.btk_version}",
                params.btk_nf_params,
            ].join(" "),
            BTK_INPUT(
                ch_reference_hap1,
                ch_btk_un_diamond_db,
                ch_btk_nt_db,
                ch_btk_un_diamond_db,
                ch_btk_ncbi_taxonomy_path,
                ch_busco_lineages,
                ch_btk_taxid,
                ch_busco_config,
                [
                    'accession': 'GCA_0001',
                    'use_work_dir_as_temp': true,
                    'align': true,
                ],
            ).json_params_file,
            GENERATE_SAMPLESHEET.out.csv,
            params.btk_extra_config ? file(params.btk_extra_config, checkIfExists: true) : [],
            workflow.workDir.resolve('sanger-tol/blobtoolkit').toUriString(),
        )
        ch_versions = ch_versions.mix(
            SANGER_TOL_BTK.out.outdir.map { outdir -> outdir.resolve('pipeline_info/sanger-tol_blobtoolkit_software_versions.yml') }
        )
    }


    //
    // LOGIC: STEP TO STOP CURATION_PRETEXT RUNNING IF SPECIFIED BY USER
    //
    if (!exclude_steps.contains("cpretext")) {

        //
        // MODULE: Run SANGER-TOL/CurationPretext
        //
        SANGER_TOL_CPRETEXT(
            'sanger-tol/curationpretext',
            [
                "-profile ${workflow.profile}",
                "-r ${params.cpretext_version}",
                params.cpretext_nf_params,
            ].join(" "),
            CPRETEXT_INPUT(
                ch_reference_hap1,
                ch_longread_dir,
                ch_cpretext_hic_dir,
                ch_cpretext_telomotif.map { it -> it[1] },
                ch_cpretext_aligner,
                [:],
            ).json_params_file,
            ch_reference_hap1,
            params.cpretext_extra_config ? file(params.cpretext_extra_config, checkIfExists: true) : [],
            workflow.workDir.resolve('sanger-tol/curationpretext').toUriString(),
        )
        ch_versions = ch_versions.mix(
            SANGER_TOL_CPRETEXT.out.outdir.map { outdir -> outdir.resolve('pipeline_info/sanger-tol_curationpretext_software_versions.yml') }
        )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ear_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
