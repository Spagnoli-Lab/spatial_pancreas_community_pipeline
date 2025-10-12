nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: "${projectDir}/samplesheet.csv"
params.outdir      = params.outdir      ?: "${projectDir}/outputs/nf"
params.sc_h5ad     = params.sc_h5ad     ?: "${projectDir}/test_data/single_cell/Cartana_simplified_fixname.h5ad"

process FILTER_CELLS_BY_TISSUE {
  tag { sample_id }
  publishDir "${params.outdir}/cell_filter", mode: 'copy'
  conda "${projectDir}/envs/filter.yml"

  input:
  tuple val(sample_id), path(mask_dir), path(cellpose_dir), val(masked_reads_path)

  output:
  tuple val(sample_id),
        path("${sample_id}_processed_mask.npz"),
        path("${sample_id}_cell_masks.tif"),
        val(masked_reads_path)

  script:
  """
  python ${projectDir}/scripts/cell_segmentation_filter.py \
    --sample ${sample_id} \
    --cellpose_path ${cellpose_dir} \
    --mask_path ${mask_dir} \
    --output_dir .
  """
}

process OVERLAY_MASKED_READS {
  tag { sample_id }
  publishDir "${params.outdir}/overlay_qc", mode: 'copy'
  conda "${projectDir}/envs/overlay.yml"

  input:
  tuple val(sample_id), path(mask_npz), path(mask_viz), path(masked_reads)

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

process RUN_PCISEQ {
  tag { sample_id }
  publishDir "${params.outdir}/pciseq", mode: 'copy'
  conda "${projectDir}/envs/pciseq.yml"

  input:
  tuple val(sample_id), path(mask_npz), path(mask_viz), path(masked_reads), path(sc_reference)

  output:
  tuple val(sample_id), path(mask_npz), path(mask_viz), path(masked_reads)
  path "${sample_id}_pciseq_result.pkl"

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

workflow {
  def sample_rows = Channel                                             // Holds rows from the samplesheet
    .fromPath(params.samplesheet)                                       // Read the CSV referenced by --samplesheet
    .splitCsv(header: true)                                             // Parse into Maps keyed by column names
    .map { row ->                                                       // Normalize each row
      def sampleId = row.sample_id                                      // Prefer an explicit sample_id column
      if (!sampleId) {                                                  // If missing, synthesize one
        def stage = (row.sample_stage ?: row.stage)                     // Pick a stage value from available columns
        def index = (row.mouse_index ?: row.sample_index ?: row.index)  // Pick an index value from available columns
        if (stage && index) {                                           // Only allow rows that can form a composite id
          sampleId = "${stage}_${index}"                                // Build a sample_id from stage/index
          row.sample_id = sampleId                                      // Store it back on the row for downstream steps
        } else {                                                        // If required pieces are absent
          throw new IllegalArgumentException(                           // Abort the run with a helpful message
            "Samplesheet row must include sample_id or stage/index columns"
          )
        }
      }
      row                                                               // Emit the normalized row into the channel
    }

  def sc_reference = file(params.sc_h5ad)                               // Load precomputed scRNA reference (AnnData)
  if (!sc_reference.exists()) {                                         // Ensure reference file is present
    throw new IllegalArgumentException("scRNA reference AnnData not found: ${sc_reference}")
  }
  log.info "[workflow] Using scRNA AnnData reference: ${sc_reference}"   // Report which reference is used

  def filtered_masks_ch = sample_rows                                    // Start from normalized sample rows
    .map { row ->                                                       // Prepare inputs for cell filtering
      def mask = file(row.mask_path)                                    // Convert mask path string to a Nextflow file
      def seg  = file(row.seg_path)                                     // Convert segmentation path likewise
      tuple(                                                            // Emit tuple expected by FILTER_CELLS_BY_TISSUE
        row.sample_id as String,
        mask.parent,
        seg.parent,
        row.reads_masked_path                                           // Keep masked reads path for downstream overlay
      )
    }
    | FILTER_CELLS_BY_TISSUE                                            // Run the cell-filtering process

  filtered_masks_ch                                                     // Prepare inputs for pciSeq (and pass-through for overlay)
    .map { sample_id, mask_npz, mask_viz, masked_reads_path ->
      tuple(
        sample_id,
        mask_npz,
        mask_viz,
        file(masked_reads_path),
        file(params.sc_h5ad)
      )
    }
    | RUN_PCISEQ

  RUN_PCISEQ.out[0]                                                     // Pass-through tuple for overlay generation
    | OVERLAY_MASKED_READS

  // Channel.value(1) | PROCESS_SCRNA                                   // Downstream scRNA step (currently disabled)
}
