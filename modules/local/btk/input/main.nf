process BTK_INPUT {

    input:
    tuple val(meta), val(reference)
    val blastp
    val blastn
    val blastx
    val tax_dump
    val busco_lineages
    val taxon
    val gca_accession
    val busco_config
    val btk_extra_opts // Map

    output:
    path "btk_params_file.json", emit: json_params_file

    exec:
    def btk_inputs = btk_extra_opts + [
        'fasta': reference,
        'busco_lineages': busco_lineages,
        'taxon': taxon,
        'taxdump': tax_dump,
        'blastp': blastp,
        'blastn': blastn,
        'blastx': blastx,
        'accession': gca_accession,
        'use_work_dir_as_temp': true,
        'align': true,
    ]
    def jsonBuilder = new groovy.json.JsonBuilder(btk_inputs)
    file("${task.workDir}/btk_params_file.json").text = jsonBuilder.toPrettyString()
}
