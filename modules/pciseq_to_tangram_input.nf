process PCISEQ_TO_TANGRAM_INPUT {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/pciseq", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(pciseq_pickle)

  output:
  tuple val(sample_id), path("${sample_id}_pciseq_predicted.h5ad")

  script:
  """
  python ${projectDir}/scripts/pciseq_to_anndata.py \\
    --sample ${sample_id} \\
    --pciseq_result ${pciseq_pickle} \\
    --output_dir .
  """
}
