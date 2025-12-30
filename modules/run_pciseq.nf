process RUN_PCISEQ {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/pciseq", mode: 'copy', saveAs: { filename ->
    def name = filename.toString()
    name.endsWith('_cell_masks.tif') ? null : filename
  }
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(mask_npz), path(mask_viz), path(masked_reads), path(sc_reference)

  output:
  tuple val(sample_id), path(mask_npz), path(mask_viz), path(masked_reads), path("${sample_id}_pciseq_result.pkl")

  script:
  """
  python ${projectDir}/scripts/run_pciseq_fit.py \\
    --sample ${sample_id} \\
    --masked_reads_csv ${masked_reads} \\
    --mask_array ${mask_npz} \\
    --sc_h5ad ${sc_reference} \\
    --cluster_key active.ident \\
    --collapse_immune \\
    --output_dir .
  """
}
