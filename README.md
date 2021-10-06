# vcfbub

popping bubbles in vg deconstruct VCFs

## overview

The VCF output produced by a command like `vg deconstruct -e -a -H '#' ...` includes information about the nesting of variants.
With `-a, --all-snarls`, we obtain not just the top level bubbles, but all nested ones.
This exposed snarl tree information can be used to filter the VCF to obtain a set of non-overlapping sites (n.b. "snarl" is a generic model of graph bubbles including tips and loops).

`vcfbub` lets us do two common operations on these VCFs:

1. We can filter sites by maximum level in the snarl tree. For instance, `--max-level 0` would keep only sites with `LV=0`. In practice, vg's snarl finder ensures that these are sites rooted on the main linear axis of the pangenome graph. Those at higher levels occur within larger variants.
2. We can filter sites by maximum allele size. In this case, `--max-length 10000` would keep only sites where the largest allele is less than 10kb. Setting `--max-length` additionally ensures that the output contains the bubbles nested inside of any popped bubble, even if they are at greater than `--max-level`.

`vcfbub` accomplishes a simple task: we keep sites that are the children of those which we "pop" due to their size.
These occur around complex large SVs, such as multi-Mbp inversions and segmental duplications.
We often need to remove these, as they provide little information for many downstream applications, such as haplotype panels or other imputation references.

## usage

This removes all non-top-level variant sites (`-l 0`) unless they are inside of variants > 10kb (`-b 10000`):

```
vcfbub -l 0 -b 10000 var.vcf >filt.vcf
```
