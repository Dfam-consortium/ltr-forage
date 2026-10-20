#!/bin/bash
# Print the 20p<GC>g scoring matrix for a FASTA file's GC content (N excluded), by RepeatMasker's
# chooseMatrices rule: the odd m with m-1 < GC <= m+1, within 35..53 (data/matrices).
# usage: pick_matrix.sh <seq.fa>
set -eo pipefail
D=$(cd "$(dirname "$0")/../../../data/matrices" && pwd)
awk -v d="$D" '!/^>/{s=toupper($0); n+=length(s); g+=gsub(/[GC]/,"",s); x+=gsub(/N/,"",s)}
END { gc = 100 * g / (n - x); h = gc / 2; c = (h == int(h)) ? h : int(h) + 1; m = 2 * c - 1
      if (m < 35) m = 35; if (m > 53) m = 53; printf "%s/20p%dg.matrix\n", d, m }' "$1"
