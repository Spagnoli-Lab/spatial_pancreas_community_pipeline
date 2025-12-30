process RUN_TANGRAM {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/tangram", mode: 'copy', saveAs: { filename ->
    def name = filename.toString()
    name.endsWith('_tangram_map.h5ad') ? null : filename
  }
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(pciseq_anndata), path(sc_reference)

  output:
  tuple val(sample_id),
        path("${sample_id}_tangram_predicted.h5ad"),
        path("${sample_id}_projected.h5ad"),
        path("${sample_id}_tangram_map.h5ad")

  script:
  """
  python ${projectDir}/scripts/run_tangram.py \\
    --sample ${sample_id} \\
    --sc_reference "${sc_reference}" \\
    --pciseq_anndata "${pciseq_anndata}" \\
    --output_dir .
  """
}
