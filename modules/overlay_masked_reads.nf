process OVERLAY_MASKED_READS {
  tag { sample_id }
  publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/overlay_qc", mode: 'copy'
  conda "${projectDir}/envs/general.yml"

  input:
  tuple val(sample_id), path(mask_npz), path(masked_reads)

  output:
  path "${sample_id}_masked_reads_overlay.png"
  path "${sample_id}_masked_reads_summary.txt"

  script:
  def sampleParts = sample_id.tokenize('_')
  def overlayStage = sampleParts ? sampleParts[0] : sample_id
  def overlayIndex = sampleParts.size() > 1 ? sampleParts[-1] : ''

  """
  python ${projectDir}/scripts/overlay_masked_reads.py \
    --sample ${sample_id} \
    --sample_stage ${overlayStage} \
    --mouse_index ${overlayIndex} \
    --masked_reads_csv ${masked_reads} \
    --mask_npz ${mask_npz} \
    --output_dir .
  """
}
