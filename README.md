# HiTea
computational tool to identify trasposable element insertions using Hi-C data


- Description:

HiTea (Hi-C based Transposable Element Analyzer) is geared to idenify non-reference transposable elemenet insertions from the Hi-C data. It uses split Hi-C read information and read coverage to detect insertions of three major classes of active human transposable elements (TE) - Alu (SINE), L1HS(LINE) and SVA. 

 - Required inputs:
HiTtea requires following reference data files during the run.  
   - TE-family consensus                = TE-family consensus sequences in fasta format. These sequences can be derived from literature or by performing multiple sequence alignment of highly variable sequences of a TE-family in the human genome. In addition, we use a PolyA fasta entry in this file to score for TPRT-mediated polyA/T tails on the non-reference TE-insertion (ref-4). HiTea uses TE-family consensus sequences for Alu, L1Hs and SVA from an earlier study (ref-5).  
   - TE-family reference copy locations = This .gzip file contains genomic location of the TE-family in bed format. The 7th column in this file is reserved for name of the TE-family as provided in the above TE-family consensus fasta file. This file can be downloded from Repeatmasker database. HiTea uses the file as provided by (ref-5) 
  
```
 $ zcat hg38/bgRepeats_hg38.bed.gz | head
chr1    26790   27053   AluSp   .       +       Alu
chr1    31435   31733   AluJo   .       +       Alu
chr1    33465   33509   Alu     .       +       Alu
chr1    35366   35499   AluJr   .       +       Alu
chr1    39623   39924   AluSx   .       +       Alu
chr1    40628   40729   AluSz6  .       -       Alu
(chr) (start) (end) (name)  (score) (strand)  (TE-family)
```
- Optional inputs:   
   - Repbase subfamily consensus      = HiTea uses Repbase subfamily consensus sequences to determine the subfamily annotation of the non-reference TE insertions
   - Polymorphic sequence file        = If the clip-sequences from Hi-C non-conforming reads do not map to the TE-family consensus sequences, the users can provide polymorphic sequences for a given TE-family. It is possible that the actual sequence of the TE-element may be divergent from the family-wise consensus. The file contains fasta format entries. The fasta header should be in the format "TE-family-name~Polymorphic-sequence-name. At present, HiTea uses Repbase (ref-6) subfamily sequences as polymorphic sequences. This file required when remap option is turned on (-remap 'T')
```
>Alu~AluYb9
ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggtggatcatgaggtc
```
 
- Running HiTea:
*It is required to set the path-variable for HiTea bash script file before running it.
  
 HiTea can be run as a single line command. First, it checks for all dependancies and once satisfied, runs on the user-provided bam file. HiTea accepts following input files
  - Lossless bam file: These name-sorted bam files can be downloaded from the 4DN data portal where the read-pairs carry classification tags generated by Pairtools [https://github.com/mirnylab/pairtools]
  - Alternately, HiTea can use name sorted bam files as an input and employs Pairtools for read-type classification
  - If the user wants to use hg19 or hg38 genome assemblies the precomputated files are available. User will need to input the genome build and RE motif information
  - At present, HiTea supports following RE enzymes: MboI/DpnII/HindIII/NcoI/Arima-cocktail
  - HiTea is currently tested on 100bp paired end sequencing data.
  
```
$ hitea
****No input provided


Usage: hitea [-w workdir] [-e enzyme] [-q anchor_mapq] [-o outprefix] [-s clip] [-g genome] [-r remap] [-n index] [-b repbase] [-p indexP] [-a anno] [-h help] -i inputs (space separated psam/bam in inverted commas)

Required****
    -i inputs :          Input file in pairsam format or unsorted-lossless bam format
    -e enzyme :          Restriction endunuclease used for the assay (default: '', available:MboI,DpnII,HindIII,Arima,NcoI,NotI)
    -g genome :          Genome build to be used (default:hg38, available: hg19)

Optional
  (following 4 parameters are optional if -g is specified)
    -n index :           fasta format file for TE-consensus sequences
    -b repbase :         fasta format file for Repbase subfamily sequences
    -p indexP :          fasta format file for Polymorphic sequences (header should be Family~name format) (optional)
    -a anno :            reference-genome copies for TE-family members

    -o outprefix :       Output prefix while generating report files (default: project)
    -w workdir:          Working directory where the files are to be written
    -q anchor_mapq :     Mapping quality threshold for repeat anchored mate on the reference genome (default: 28)
    -s clip :            Minimum clip length for detecting insertion (should be >=13bp) (default: 20)
    -r remap :           whether to remap unmapped clipped reads to the polymoprhic sequences (default:F, if T, -p needs to be specified)
    -h help :            Display help message


```
  
Running HiTEA on a single input bam file
```
hitea -i bam/4DNFIC275NK8.bam -w GM12878_test -o gm12878 -g hg38 -e 'MboI' -r 'T'
```

Running HiTEA on a list of input bam files of a single experiment
```
hitea -i 'bam/4DNFPC275NK8.bam bam/4DNFIJ275PQ9.bam bam/4DNFIC275HT2.bam' -w GM12878_test -o gm12878 -g hg38 -e 'MboI' -r 'T'
```

# Installation
- Dependancies:
  - PERL(≥v5.24)
  - R(≥v3.2)
  - bedtools(≥v2.26) (ref-1)
  - samtools(≥v1.7)
  - GNU-parallel [https://www.gnu.org/software/parallel/] (ref-2)
  - Pairtools [https://github.com/mirnylab/pairtools]
  - Following R packages are required for computation
    - GenomicRanges, data.table, MASS
  - Following optional R packages are required for HTML report generation
    - rmarkdown, knitr, EnrichedHeatmap (ref-3), circlize
 
 - HiTea can be installed in following three ways:
 
- 1) By setting up the environment:
    - Once the dependancies are installed, the HiTea package can be installed simply by downloading it and setting the path variable to the HiTea directory. 

- 2) Docker image: 
    - HiTea docker image is available through 4DN DCIC DockerHub repo (4dndcic/hitea:v1)

- 3) Conda package:
    - Additionally, HiTea is available through bioconda. 
    - Highly recommended option to run HiTea on local machines
```
conda install -c bioconda hitea
```


# Codes used for analyses
- Codes used for analyzing the data and preparing figures for manuscript are located here: 
https://github.com/dhawalsjain/HiTea_Analyses


References:
1. Quinlan,A.R. and Hall,I.M. (2010) BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26, 841–842
2. Tange,O. (2011) GNU Parallel: The Command-Line Power Tool. USENIX Mag., 36, 42–47
3. Gu,Z. et al. (2018) EnrichedHeatmap: An R/Bioconductor package for comprehensive visualization of genomic signal associations. BMC Genomics, 19, 234
4. Lee,E. et al. (2012) Landscape of somatic retrotransposition in human cancers. Science (80-.). doi:10.1038/nrg2072
5. Gardner,E.J. et al. (2017) The mobile element locator tool (MELT): Population-scale mobile element discovery and biology. Genome Res., 27, 1916–1929. doi:10.1101/gr.218032.116
6. Bao,W. et al. (2015) Repbase Update, a database of repetitive elements in eukaryotic genomes. Mob. DNA, 6:11
