nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: "${projectDir}/samplesheet.csv"
params.outdir      = params.outdir      ?: "${projectDir}/outputs/nf"
params.stage       = params.stage       ?: "E14"
params.sc_expr     = params.sc_expr     ?: "${projectDir}/test_data/single_cell/CartanaJF_singlecell_df.txt"
params.sc_meta     = params.sc_meta     ?: "${projectDir}/test_data/single_cell/CartanaJF_singlecell_meta_df.txt"

process FILTER_CELLS_BY_TISSUE {
  tag { sample_id }
  publishDir "${params.outdir}/cell_filter", mode: 'copy'
  conda "${projectDir}/envs/filter.yml"

  input:
  tuple val(sample_id), path(mask_dir), path(cellpose_dir)

  output:
  path "${sample_id}_processed_mask.npz"
  path "${sample_id}_cell_masks.tif"

  script:
  """
  python ${projectDir}/scripts/cell_segmentation_filter.py \
    --sample ${sample_id} \
    --cellpose_path ${cellpose_dir} \
    --mask_path ${mask_dir} \
    --output_dir .
  """
}

process PROCESS_SCRNA {
  publishDir "${params.outdir}/scRNA", mode: 'copy'
  conda "${projectDir}/envs/scrna.yml"

  input:
  val x

  output:
  path "sc_${params.stage}.h5ad"

  script:
  """
  python ${projectDir}/scripts/process_scdata.py \
    --scRNA_df ${params.sc_expr} \
    --scRNA_meta_df ${params.sc_meta} \
    --stage ${params.stage} \
    --output_h5ad sc_${params.stage}.h5ad
  """
}

workflow {
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      def mask = file(row.mask_path)
      def seg  = file(row.seg_path)
      tuple(row.sample_id as String, mask.parent, seg.parent)
    }
    | FILTER_CELLS_BY_TISSUE

  Channel.value(1) | PROCESS_SCRNA
}
