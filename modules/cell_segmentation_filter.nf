process CELL_SEGMENTATION_FILTER {
    tag { sample_id }
    publishDir "${params.outdir ?: projectDir + '/outputs/nf'}/${sample_id}/cell_filter", mode: 'copy'
    conda "${projectDir}/envs/general.yml"

    input:
    tuple val(sample_id), path(cellpose_dir), path(mask_dir), val(masked_reads_path)

    output:
    tuple val(sample_id),
          path("${sample_id}_processed_mask.npz"),
          path("${sample_id}_cell_masks.tif"),
          val(masked_reads_path)

    script:
    """
    python ${projectDir}/scripts/cell_segmentation_filter.py \\
        --sample ${sample_id} \\
        --cellpose_path ${cellpose_dir} \\
        --mask_path ${mask_dir} \\
        --output_dir .
    """
}
