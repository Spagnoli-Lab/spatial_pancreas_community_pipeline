# Nextflow pipeline arguments — explanation and examples

This document explains the main Nextflow runtime options and pipeline parameters used by the spatial_pancreas_community_pipeline, plus recommended defaults and example invocations. It is intended to help maintainers and users confirm argument settings and to serve as the canonical reference for reproducible runs.

## Location

File: docs/nextflow-arguments.md

## Purpose

- Confirm which pipeline arguments are expected by the workflow.
- Explain default values (where applicable) and recommended settings for common environments (local, HPC/Slurm, cloud).
- Provide example commands and a sample params YAML for common use-cases.

## Nextflow runtime options (commonly used)

- -profile <profile>
  - Selects a Nextflow profile defined in nextflow.config (e.g., 'standard', 'docker', 'singularity', 'slurm', 'aws').
  - Recommendation: Use explicit profile names rather than relying on the implicit defaults.

- -resume
  - Resume a previous run by reusing cached work. Use when re-running after fixing parameters or inputs.

- -with-dag <file>.dot / -with-trace / -with-timeline
  - Useful for debugging and provenance.

- -c <config-file>
  - Use a custom Nextflow config file if you need site-specific executor settings.

## Pipeline parameters (params.*)

Below are the parameters the pipeline typically accepts. Replace names with the exact parameter keys used by the workflow (check main.nf or nextflow.config if uncertain).

- params.input_dir (string)
  - Path to directory containing input data (e.g., fastq/, images/). Recommended: absolute paths or use S3/GS URI when running in cloud.

- params.sample_sheet (string)
  - CSV/TSV describing samples and metadata. If provided, the pipeline will use this to discover inputs rather than scanning input_dir.

- params.genome (string)
  - Reference genome identifier or path (e.g., mm10, hg38, or /path/to/reference.fasta).

- params.gtf (string)
  - Path to annotation GTF file when needed by the pipeline.

- params.outdir (string)
  - Output directory for pipeline results. Default: results/ or workspace/; recommendation: provide an explicit path.

- params.max_cpus (int) and params.max_memory (string)
  - Upper resource bounds for the pipeline; these act as safety limits.

- params.threads (int)
  - Default number of threads per task if not overridden in process config.

- params.container_image (string)
  - If using Docker/Singularity, the container image name (and tag) to use. Pin to a digest or exact tag for reproducibility.

- params.skip_* (boolean)
  - Flags to skip optional pipeline stages (e.g., skip_qc, skip_alignment).

- params.extra (map)
  - A flexible map for site-specific overrides.

Note: The exact parameter names may vary. Open the pipeline's main.nf and nextflow.config to verify names; this document should be updated if parameter names change.

## Example params.yaml

```yaml
input_dir: '/data/project/inputs'
sample_sheet: '/data/project/samples.csv'
genome: 'hg38'
gtf: '/data/reference/hg38.gtf'
outdir: '/data/project/results'
max_cpus: 64
max_memory: '256.GB'
container_image: 'ghcr.io/spagnoli-lab/spatial-pancreas:2025-10-01'
skip_qc: false
```

## Example commands

- Local run using Docker (for reproducible local testing):

```bash
nextflow run . \
  -profile docker \
  -params-file params.yaml \
  -resume
```

- HPC run using Singularity + Slurm profile:

```bash
nextflow run . \
  -profile slurm,singularity \
  -params-file params.yaml
```

- Cloud run (AWS Batch example):

```bash
nextflow run . \
  -profile awsbatch \
  -params-file params.yaml
```

## Recommended checks to confirm argument settings

1. Verify parameter names in main.nf and nextflow.config match what users expect (search for "params.").
2. Confirm default values in nextflow.config and document them here.
3. Confirm profiles in nextflow.config (names and what they set: executor, container engines, queue names, memory/cpu defaults).
4. Confirm container image tags and whether images are published to a registry; pin images for reproducibility.
5. Add a CI job or a checklist in the repository to verify the above when merging changes that alter defaults or profiles.

## Troubleshooting tips

- If tasks fail due to memory or CPU, increase params.max_memory or params.max_cpus, and consider adjusting process-level resource declarations.
- For permission or runtime errors when using Singularity on an HPC, ensure the singularity cache directory is writable by the user and set SINGULARITY_CACHEDIR if needed.
- If reproducibility is required, avoid using :latest tags for container images.

## Next steps

- Maintainers: please confirm the exact parameter names and defaults by reviewing main.nf and nextflow.config and update this file accordingly.
- Once confirmed, this file should be kept in docs/ and referenced from the README.

---

(End of file)