// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

process GTF_REMOVE_DIRT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.clean.gtf")    , emit: gtf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'NR==FNR {
            if (\$4 >= \$5) { split(\$10, a, "\\""); genes[a[2]] = 1 }
            next
        }
        {
            split(\$10, a, "\\"")
            if (!(a[2] in genes)) print
        }' ${gtf} ${gtf} > ${prefix}.clean.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | sed 's/^.*gawk version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clean.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | sed 's/^.*gawk version //; s/ .*//')
    END_VERSIONS
    """
}               
