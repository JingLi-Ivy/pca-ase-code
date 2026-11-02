#!/usr/bin/env bash
# Point the placeholder paths in the scripts to your own data layout.
#
# The scripts refer to inputs and outputs under the placeholder root /path/to/data_root
# (this stands for the directory layout used when the reported results were produced).
# This helper rewrites the placeholder on a COPY of the repository.
#
# Usage:
#   cp -r pca-ase-code pca-ase-code-local
#   ./config/remap_paths.sh pca-ase-code-local /your/data/root
#
# Expected layout under /your/data/root:
#   ASE/            <- project root
#   CPGEA/          <- expression matrices and clinical table
#   reference/      <- genome FASTA and annotation
#   raw/            <- controlled-access sequencing data (read-only inputs)
#   variants/       <- Mutect2 and FACETS outputs
#   singlecell/     <- single-cell reference and deconvolution output
#   wasp/           <- STAR/WASP working directory
#   user/           <- other per-user inputs (e.g. the external ASE caller checkout)
# /path/to/software/share is the location of shared tools (e.g. git modules); edit by hand.

set -euo pipefail

TARGET="${1:?usage: remap_paths.sh <repo-copy-dir> <new-root>}"
ROOT="${2:?usage: remap_paths.sh <repo-copy-dir> <new-root>}"

[ -d "$TARGET" ] || { echo "not a directory: $TARGET" >&2; exit 1; }

find "$TARGET" -type f \( -name '*.R' -o -name '*.sh' -o -name '*.py' \) -print0 |
while IFS= read -r -d '' f; do
  ROOT="$ROOT" perl -pi -e 's{/path/to/data_root}{$ENV{ROOT}}g' "$f"
done

echo "Remapped /path/to/data_root to $ROOT under $TARGET"
echo "Remaining placeholders (review manually):"
grep -rn "/path/to/" "$TARGET" --include='*.R' --include='*.sh' --include='*.py' | head -20 || true
