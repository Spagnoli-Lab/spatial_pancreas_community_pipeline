#!/usr/bin/env python3
"""Generate masked DAPI images by combining 8-bit masks with their raw DAPI TIFFs.

The utility defaults to the repository test data layout where masks live under
``test_data/8_bit_masks`` and the matching DAPI images under
``test_data/8_bit_masks/dapi``. For each match the script now emits two outputs:

* ``full`` – masked image with the original spatial dimensions preserved
* ``cut_black_edge`` – masked image cropped to the minimal bounding box around
  non-zero pixels (removing the black background)

Use the CLI flags to override locations when working with other datasets.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np

try:
    import tifffile
except ImportError as exc:  # pragma: no cover - ensure a helpful message if dependency is missing
    raise SystemExit("The 'tifffile' package is required. Install it with `pip install tifffile`.") from exc


MASK_TOKEN = re.compile(r"(?i)(dapi_)?(e\d+(?:\.\d+)?)_(\d+)")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mask-dir",
        type=Path,
        default=Path("test_data") / "8_bit_masks",
        help="Directory containing 8-bit mask TIFFs (default: test_data/8_bit_masks)",
    )
    parser.add_argument(
        "--dapi-dir",
        type=Path,
        default=Path("test_data") / "8_bit_masks" / "dapi",
        help="Directory containing raw DAPI TIFFs (default: test_data/8_bit_masks/dapi)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("test_data") / "8_bit_masks" / "masked_dapi",
        help="Destination for masked DAPI TIFFs (default: test_data/8_bit_masks/masked_dapi)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs when they already exist",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Inspect planned work without writing output files",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args(argv)


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(message)s")


def iter_mask_files(mask_dir: Path) -> Iterable[Path]:
    if not mask_dir.exists():
        raise SystemExit(f"Mask directory does not exist: {mask_dir}")
    yield from sorted(p for p in mask_dir.glob("*.tif") if p.is_file())


def sample_id_from_mask(path: Path) -> Optional[str]:
    match = MASK_TOKEN.search(path.stem)
    if not match:
        return None
    stage = match.group(2)
    index = match.group(3)
    return f"{stage.upper()}_{index}"


def find_matching_dapi(dapi_dir: Path, sample_id: str) -> Path:
    if not dapi_dir.exists():
        raise SystemExit(f"DAPI directory does not exist: {dapi_dir}")

    candidates = [
        dapi_dir / f"DAPI_{sample_id}.tif",
        dapi_dir / f"{sample_id}.tif",
        dapi_dir / f"{sample_id.lower()}.tif",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    available = sorted(p.name for p in dapi_dir.glob("*.tif"))
    raise SystemExit(
        f"Could not find a DAPI TIFF for sample '{sample_id}' in {dapi_dir}. "
        f"Checked candidates: {', '.join(c.name for c in candidates)}. "
        f"Available DAPI files: {', '.join(available) if available else 'none'}"
    )


def load_mask(mask_path: Path) -> np.ndarray:
    logging.debug("Loading mask %s", mask_path)
    return tifffile.imread(mask_path)


def load_dapi(dapi_path: Path) -> np.ndarray:
    logging.debug("Loading DAPI %s", dapi_path)
    return tifffile.imread(dapi_path)


def normalize_mask(mask: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    mask_bool = np.squeeze(mask > 0)
    if mask_bool.shape == target_shape:
        return mask_bool

    if mask_bool.shape == target_shape[:-1]:
        mask_bool = mask_bool[..., np.newaxis]

    try:
        mask_bool = np.broadcast_to(mask_bool, target_shape)
    except ValueError as exc:
        raise SystemExit(
            f"Mask shape {mask_bool.shape} is incompatible with DAPI shape {target_shape}"
        ) from exc
    return mask_bool.astype(bool, copy=False)


def crop_to_mask(image: np.ndarray, mask_bool: np.ndarray) -> np.ndarray:
    if not np.any(mask_bool):
        return image

    reduced = mask_bool
    while reduced.ndim > 2:
        reduced = np.any(reduced, axis=-1)

    coords = np.argwhere(reduced)
    if coords.size == 0:
        return image

    mins = coords.min(axis=0)
    maxs = coords.max(axis=0) + 1

    spatial_slices = [slice(start, stop) for start, stop in zip(mins, maxs)]
    extras = [slice(None)] * (image.ndim - len(spatial_slices))
    cropped = image[tuple(spatial_slices + extras)]
    return np.ascontiguousarray(cropped)


def apply_mask_to_dapi(dapi: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask_bool = normalize_mask(mask, dapi.shape)
    masked_full = np.zeros_like(dapi)
    masked_full[mask_bool] = dapi[mask_bool]
    masked_cut = crop_to_mask(masked_full, mask_bool)
    return masked_full, masked_cut


def process_sample(
    sample_id: str,
    mask_path: Path,
    dapi_path: Path,
    output_dir: Path,
    overwrite: bool,
    dry_run: bool,
) -> None:
    full_dir = output_dir / "full"
    cut_dir = output_dir / "cut_black_edge"
    full_dir.mkdir(parents=True, exist_ok=True)
    cut_dir.mkdir(parents=True, exist_ok=True)

    full_output = full_dir / f"{sample_id}_masked_dapi.tif"
    cut_output = cut_dir / f"{sample_id}_masked_dapi.tif"

    if not overwrite and full_output.exists() and cut_output.exists():
        logging.info("Skipping %s (exists); use --overwrite to regenerate", sample_id)
        return

    if dry_run:
        logging.info(
            "[DRY-RUN] Would create %s and %s using %s and %s",
            full_output,
            cut_output,
            dapi_path,
            mask_path,
        )
        return

    dapi = load_dapi(dapi_path)
    mask = load_mask(mask_path)
    masked_full, masked_cut = apply_mask_to_dapi(dapi, mask)

    logging.info("Writing %s", full_output)
    tifffile.imwrite(full_output, masked_full)
    logging.info("Writing %s", cut_output)
    tifffile.imwrite(cut_output, masked_cut)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    configure_logging(args.verbose)

    for mask_path in iter_mask_files(args.mask_dir):
        sample_id = sample_id_from_mask(mask_path)
        if not sample_id:
            logging.warning("Skipping %s (could not determine sample id)", mask_path.name)
            continue

        try:
            dapi_path = find_matching_dapi(args.dapi_dir, sample_id)
        except SystemExit as exc:
            logging.warning(str(exc))
            continue

        process_sample(
            sample_id=sample_id,
            mask_path=mask_path,
            dapi_path=dapi_path,
            output_dir=args.output_dir,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )


if __name__ == "__main__":
    main()
