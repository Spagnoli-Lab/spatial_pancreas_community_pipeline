process PROCESS_SCDATA {
  tag { stage }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/sc_reference", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple path(sc_expr), path(sc_meta), val(stage), val(output_name)

  output:
  path "${output_name}"

  script:
  """
  python ${projectDir}/scripts/process_scdata.py \\
    --scRNA_df ${sc_expr} \\
    --scRNA_meta_df ${sc_meta} \\
    --stage ${stage} \\
    --output_h5ad ${output_name}
  """
}
