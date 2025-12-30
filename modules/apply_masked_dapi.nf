process APPLY_MASKED_DAPI {
  tag { mask_dir }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/masked_dapi", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple path(mask_dir), path(dapi_dir)

  output:
  path("full/*_masked_dapi.tif")
  path("cut_black_edge/*_masked_dapi.tif")

  script:
  """
  python ${projectDir}/scripts/apply_masked_dapi.py \\
    --mask-dir ${mask_dir} \\
    --dapi-dir ${dapi_dir} \\
    --output-dir . \\
    --overwrite
  """
}
