#!/usr/bin/env nextflow

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METASSEMBLE } from './subworkflows/metassembly/main'
include { EMAIL_RESULTS } from './modules/custom/email/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    use_mailx
    ch_samplesheet

    main:

    if (params.sent_email) {
        EMAIL_RESULTS (
            email,
            email_on_fail,
            plaintext_email,
            outdir,
            use_mailx,
            ch_samplesheet
        )
    }

    workflow.onError {
        log.error "ERROR: Pipeline failed. Please refer to github issues: https://github.com/alejandrogzi/metassembly/issues"
    }

    workflow.onComplete {
        log.info "\nPipeline completed successfully!"
    }
}

workflow XASM {
    METASSEMBLE (
        params.input_dir,
        params.genome,
        params.annotation,
        params.star_index_path,
        params.star_ignore_gtf_for_index,
        params.deacon_index_path,
        params.deacon_download_index,
        params.deacon_make_single_index,
        params.star_make_coverage,
        params.skip_assembly,
        params.isotools_orphan_splicing_scores_dir,
        params.output_dir
    )

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.output_dir,
        params.use_mailx,
        METASSEMBLE.out.samplesheet
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    XASM()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
