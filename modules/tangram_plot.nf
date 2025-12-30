process TANGRAM_PLOT {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/tangram/post", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(spatial_adata), path(projected_h5ad), path(tangram_map), path(sc_reference)

  output:
  path "${sample_id}_tangram_training_scores.png"
  path "${sample_id}_tangram_scatter.png"
  path "${sample_id}_tangram_celltype_counts.txt"

  script:
  """
  python ${projectDir}/scripts/tangram_postprocess.py \\
    --sample ${sample_id} \\
    --spatial_adata "${spatial_adata}" \\
    --tangram_map "${tangram_map}" \\
    --output_dir .
  """
}
