nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: "${projectDir}/samplesheet.csv"
params.outdir      = params.outdir      ?: "${projectDir}/outputs/nf"
params.sc_h5ad     = params.sc_h5ad     ?: "${projectDir}/test_data/single_cell/Cartana_simplified_fixname.h5ad"

def sampleSubdir(String sampleId) {
  if (!sampleId) {
    return "unknown_sample"
  }
  def parts = sampleId.tokenize('_')
  if (parts.size() >= 2) {
    return "${parts[0]}/${parts[-1]}"
  }
  return sampleId
}

include { CELL_SEGMENTATION_FILTER } from './modules/cell_segmentation_filter.nf'
include { OVERLAY_MASKED_READS     } from './modules/overlay_masked_reads.nf'
include { RUN_PCISEQ              } from './modules/run_pciseq.nf'
include { PCISEQ_TO_TANGRAM_INPUT } from './modules/pciseq_to_tangram_input.nf'
include { RUN_TANGRAM             } from './modules/run_tangram.nf'
include { TANGRAM_PLOT            } from './modules/tangram_plot.nf'

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

  def filter_inputs = sample_rows                                       // Prepare inputs for cell filtering
    .map { row ->                                                       //
      def mask = file(row.mask_path)                                    // Convert mask path string to a Nextflow file
      def seg  = file(row.seg_path)                                     // Convert segmentation path likewise
      tuple(                                                            // Emit tuple expected by CELL_SEGMENTATION_FILTER
        row.sample_id as String,
        seg.parent,
        mask.parent,
        file(row.reads_masked_path)                                     // Carry masked reads path for downstream steps
      )
    }

  def filtered_with_reads = filter_inputs                               // Run the cell-filtering module (masked reads passthrough)
    | CELL_SEGMENTATION_FILTER

  filtered_with_reads                                                   // Generate overlay QC
    .map { sample_id, mask_npz, mask_viz, masked_reads ->
      tuple(sample_id, mask_npz, masked_reads)
    }
    | OVERLAY_MASKED_READS

  def pciseq_results = filtered_with_reads                              // Run pciSeq fit (single output tuple)
    .map { sample_id, mask_npz, mask_viz, masked_reads ->
      tuple(sample_id, mask_npz, mask_viz, masked_reads, sc_reference)
    }
    | RUN_PCISEQ

  def tangram_input = pciseq_results                                    // Convert pciSeq results into AnnData
    .map { sample_id, mask_npz, mask_viz, masked_reads, pciseq_pickle ->
      tuple(sample_id, pciseq_pickle)
    }
    | PCISEQ_TO_TANGRAM_INPUT

  def tangram_results = tangram_input                                   // Run Tangram using single-cell reference
    .map { sample_id, pciseq_anndata ->
      tuple(sample_id, pciseq_anndata, sc_reference)
    }
    | RUN_TANGRAM

  TANGRAM_PLOT(                                                         // Create Tangram visualization
    tangram_results.map { sample_id, tangram_predicted, projected_h5ad, tangram_map ->
      tuple(sample_id, tangram_predicted, projected_h5ad, tangram_map, sc_reference)
    }
  )

  // Channel.value(1) | PROCESS_SCRNA                                   // Downstream scRNA step (currently disabled)
}
