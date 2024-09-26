# Manual genome assembly purging

Procedure to manually purge genome assemblies using kmer content. Assemblers have sometimes problems to produce balanced haplotypes. When this is the case often one (or several) haplotypic assembly(ies) harbors too many contigs and the other(s) too few. It is possible using a kat kmer spectra-cn plot to visualize this state with multiple kmers in the assembly when only one copy is expected. Usually this is resolved by automatically purging the assembly within the assembler (such as in hifiasm) or with another software packag such as purge_dups. This procedure often gets rid of duplicated contigs but not always and at the moment purge_dups is not able to automatically purge contigs in triploid, tetraploid, plus assemblies. 

This manual procedure is ment to enable assembly purging before scaffolding using the read dans haplotypic assemblies kmer content. 

## tools used in this procedure 

To follow this procedure you will need
- jellyfish including query_per_sequence (https://github.com/gmarcais/Jellyfish)
- KMC (https://github.com/refresh-bio/KMC)
- kat (https://github.com/TGAC/KAT)
- seq (https://seq-lang.org/)
- python (https://www.python.org/)

## Defining the kmer groups from the kmer histogram 

In this step we will define the kmer group boundaries from the kmer histogram. In the kmer histogram we are particularly intersted in the homozygous kmer gaussian because these kmers should be present once and only once in each haplotypic assemblies. Defining homozygous kmer profile boundaries can be quite trivial if the kmer profile is made with enough coverage enabling an clean distinction of the differents overlapping gaussians or can be tricky if the coverage is limited and the gaussians are mingled. 

Here is an example of a tetraploid genome for which the kmer gaussians are well a part enable a quite simple decision about the boundaries. We fixed the thresholds from 200 to 250 coverage (blue bars) .

## Building KMC kmer databases for each read coverage 


## 
