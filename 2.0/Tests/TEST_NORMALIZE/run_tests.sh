#!/bin/bash

set -exo pipefail

# Regression test for --normalize (left-alignment) of indels that sit at the
# right edge of a tandem repeat and therefore have to be shifted left by more
# than one base.
#
# A former bug in VNormalizeContig() assumed that any allele whose original
# bases were entirely removed by right-trimming collapsed to a single reference
# base.  That's only true when the left-shift is short; when the shift spans
# multiple bases the normalized allele is several reference bases long.  The bug
# dropped everything but one base, turning insertions and deletions into
# spurious SNPs (e.g. G:GC -> G:C, GC:G -> C:G), which then silently passed
# downstream since they still looked like valid variants.
#
# normalize_test.fa is  TGCCCCAAAAAGTACGT  (1-based):
#   2=G, 3-6=CCCC (poly-C preceded by G), 6=C, 7-11=AAAAA (poly-A preceded by C)
# so all three input variants left-align several bases upstream:
#   del1C  1:5 CC>C    -> 1:2 GC>G   (deletion, longer allele is REF)
#   ins1C  1:6 C>CC    -> 1:2 G>GC   (1-base insertion)
#   ins2A  1:11 A>AAA  -> 1:6 C>CAA  (multi-base insertion)

$1/plink2 $2 $3 --vcf normalize_test.vcf --fa normalize_test.fa --normalize --make-pgen --out plink2_test
diff -q plink2_test.pvar normalize_test.want.pvar
