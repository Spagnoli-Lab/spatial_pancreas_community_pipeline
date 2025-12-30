process PLOT_TANGRAM_ANNOTATION {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/tangram/post", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(sc_reference), path(tangram_anndata)

  output:
  path "${sample_id}_tangram_annotation.png"

  script:
  """
  python ${projectDir}/scripts/plot_tangram_annotation.py \\
    --sample ${sample_id} \\
    --sc_reference ${sc_reference} \\
    --tangram_anndata ${tangram_anndata} \\
    --output_dir . \\
    --cluster_key ${params.tangram_cluster_key ?: 'active.ident'}
  """
}
