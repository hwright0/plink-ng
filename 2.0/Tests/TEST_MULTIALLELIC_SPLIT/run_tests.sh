#!/bin/bash

set -exo pipefail

# Regression test for '--make-pgen multiallelics=-' splitting of multiallelic
# variants that carry hardcall phase.  Three former bugs, all in the split
# encoder in MakePgenRobust(), are exercised here:
#
#   1. When a split output used the "explicit phasepresent" branch (a variant
#      with a mix of phased and unphased hets, so phasepresent_ct < het_ct), the
#      write pointer was not advanced past the phase track, so the next allele's
#      record overlapped it.  The reader then read a garbage het_ct and the
#      pgen writer crashed (heap-buffer-overflow / AppendHphase assertion).
#   2. The "explicit phasepresent" flag bit was set and then immediately cleared
#      by a following memset, so the reader treated the record as fully phased
#      and marked genuinely-unphased hets (e.g. S1 0/2) as phased.
#   3. Each allele's "regular" sample-index list is filled in two passes per
#      word (0/x hets first, then x/y patch_10 hets), so a patch_10 het at a
#      lower sample index than a preceding 0/x het left the list out of order.
#      The encoder walked it in list order while the reader reconstructs phase
#      in sample order, silently swapping phase between samples (e.g. S0 1|2 and
#      S2 2|0 exchanged phase direction).
#
# split_test.vcf covers, on multiallelic sites with mixed phasing:
#   - phased ref/alt and alt/alt hets in both directions (0|1, 1|2, 2|0, ...)
#   - unphased hets that must stay unphased (0/2, 1/2, 0/1)
#   - the bug-3 trigger: a low-index x/y patch_10 het plus higher-index 0/x
#     hets on the same output allele (S0 vs S1/S2 for allele G)
#   - homozygous ALT (incl. rare-allele x/x), homozygous REF, and missing calls
#
# split_test.want.txt is the hand-verified correct split: every sample's
# genotype reconstructs its original alleles, phased hets keep their haplotype
# direction, and unphased hets stay unphased.
$1/plink2 $2 $3 --vcf split_test.vcf --double-id --make-pgen multiallelics=- --out plink2_split
$1/plink2 $2 $3 --pfile plink2_split --export vcf --out plink2_split

# Drop the header lines (they include a volatile fileDate); compare the variant
# records, which carry the split genotypes and phase.
grep -v '^#' plink2_split.vcf > plink2_split_data.txt
diff -q plink2_split_data.txt split_test.want.txt
