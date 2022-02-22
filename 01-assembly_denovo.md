---
layout: page
title: Introduction to RNA-seq and functional interpretation
subtitle: De novo transcriptome assembly
minutes: 5
---

> ## Learning Objectives {.objectives}
>
> *  Learn how the *de novo* RNA-Seq assembly works.
> *  Learn how to use Trinity to assemble a transcriptome *de novo*.

To start, we will create a directory for our data:

~~~ {.bash}
mkdir De_Novo_Assembly
cd De_Novo_Assembly
~~~ 

And we will download the fastq files using wget:

~~~ {.bash}
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_ds.left.fq.gz 
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_ds.right.fq.gz
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_hs.right.fq.gz 
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_hs.left.fq.gz   
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_plat.left.fq.gz
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_plat.right.fq.gz
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_log.left.fq.gz  
$ wget https://liz-fernandez.github.io/EBI-Intro-RNA-seq-2022/datasets/Sp_log.right.fq.gz
~~~

Let's explore the reads in each file using a 'for' cycle:

~~~ {.bash}
$ for fastq in S*fq.gz ; do echo $fastq; zcat $fastq | head ; wait ; done 
~~~

This command will show us the first 10 lines of each file
in an iterative way. This structure is known as a `for loop`.
A for loop allows us to execute a command in several files
sequentially. It is extremely useful when we have several
files.

Once we verify that the data is correct,
We will use the Trinity program to assemble the transcripts.
In this exercise we are assuming that these files were already made
filtered by quality.

We will first start up conda:

~~~ {.bash}
$ conda activate
~~~

And use Trinity with a generic command:

~~~ {.bash}
$ Trinity --seqType fq --SS_lib_type RF \
--left Sp_log.left.fq.gz,Sp_hs.left.fq.gz \
--right Sp_log.right.fq.gz,Sp_hs.right.fq.gz \
--CPU 2 --no_bowtie --max_memory 2G
~~~

This command will take approximately 15 minutes to assemble a
transcriptome.

The options (flags) that we have used are the following:

* --seqType fq - We indicate that we are using fastq type files
* --SS_lib_type RF - We indicate that the library was built using
paired-end reads in the RF orientation (reverse-forward)
* --left - Left side reads (o R1)
* --right - Right side reads (o R2)
* --CPU 2 - Use 2 CPUs
* --no_bowtie  - Do not use bowtie (normally not used, but needed in this case due to the conda version of Trinity we are using)
* --max_memory 2G - Indicate that the maximum RAM to use is 2GB

These are the most essential options to carry out the analysis.
Therefore it is very important to know the orientation of the library
which depends on the experimental protocol that was used to generate it.

> ## Library strand orientation {.callout}
>
> The library orientation is an essential parameter for the
> assembly, since this will dictate the orientation of the transcripts.
>
> This is a key that will allow you to decide which parameter to use
> to indicate to Trinity the orientation of its library.
> 
> ![Library orientation](fig/strand_specificity.jpg)
>

> ## Why do we mix libraries? {.challenge}
>
> In the command used at the top we can see that
> We are assembling the transcriptome with the log library
> as well as the library hs. Why not assemble them
> independently? What advantages or disadvantages do you think it would have
> assemble them together or separately?

One of the common problems with Trinity is the lack of memory.
A typical Trinity analysis (without digital normalization) 
requires ~1 hour and ~1GB of RAM
for ~1 million paired-end reads. This is why 
these analyzes are carried out on a server with large
amounts of memory where it can be left running for several days

Trinity has many other options which we can
explore writing (in a different terminal than the one we are
using for our analysis):

~~~ {.bash}
$ Trinity
~~~

~~~ {.output}
###############################################################################
#
#     ______  ____   ____  ____   ____  ______  __ __
#    |      ||    \ |    ||    \ |    ||      ||  |  |
#    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
#    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
#      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
#      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
#      |__|  |__|\_||____||__|__||____|  |__|  |____/
#
###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( fa, or fq )
#
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provied in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
#
####################################
##  Misc:  #########################
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#
#  --CPU <int>                     :number of CPUs to use, default: 2
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=200)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format
#                                   for reads).
#                                   (note: jaccard_clip is an expensive
#                                   operation, so avoid using it unless
#                                   necessary due to finding excessive fusion
#                                   transcripts w/o it.)
#
#  --trimmomatic                   :run Trimmomatic to quality trim reads
#                                        see '--quality_trimming_params' under full usage info for tailored settings.
#                                  
#
#  --normalize_reads               :run in silico normalization of reads. Defaults to max. read coverage of 50.
#                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
#     
#  --no_distributed_trinity_exec   :do not run Trinity phase 2 (assembly of partitioned reads), and stop after generating command list.
#
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/home/training/Documents/De_Novo_Assembly/trinity_out_dir" 
#                                    note: must include 'trinity' in the name as a safety precaution! )
#  
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (v2.1.1) and exits.
#
#  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
#
#
###############################################################################
#
#  *Note, a typical Trinity command might be:
#
#        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
#
#
#    and for Genome-guided Trinity:
#
#        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
#
#     see: /home/training/anaconda3/opt/trinity-2.1.1/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#
#     For more details, visit: http://trinityrnaseq.github.io
#
###############################################################################
~~~ 

Our work must be finished, let's review the transcripts
generated, which are in the file `Trinity.fasta`:

~~~ {.bash}
$ head trinity_out_dir/Trinity.fasta
~~~

We note that the results are sequences in `fasta` format.

~~~ {.output}
>TRINITY_DN8_c0_g1_i1 len=1720 path=[1:0-1719] [-1, 1, -2]
TTGCAATGCAAGTATTTAAGCGTATCACAACACATTGTTCTTCTCCAAGGTTCGTGAATC
GCTGTATTATTCTTTATTTTTCTTCAAGTGAGGATAAGAGTGATTGTTTGGCAAAGAAAA
ATTATGTTAACAAGTGTCTTATGGCCAAGGCACTTAAAGATTACCCCGTTCATACCAATA
TTGATCCTGATGCAGGGAAGTTATCATTTGACGATGCTTTTTACGAAGCTCACATTGAAC
TTCATTATCAATTTTTGAAGGAGGCTTCCCTAAATACCCTTATTAAAGATAAAAAAATGC
TCAAGTTTATTATCACTGTTCGTCCCGTTCATTTGCATGTCTCACCTTGGGTGGTTTATC
GTCGATATCGGGGGTTCAAAACTTTATACTATTTGTTAAAAAAGCAAAGTGCTAGAAATG
GGCGAGCTGTACCGAGTTTTCCTGTTTGGCGTGGAAACACGTATGAGAAGTTTCGTGAAG
GATTGTATTTTTTTATAGAAGCTTTACTGCATGATAGTCACTTTGCAACTAATGTTGATG
~~~

## Analyzing the statistics of the assembled transcriptome

We can capture some statistics about this assembly
using a program that is part of Trinity:

~~~ {.bash}
$ TrinityStats.pl trinity_out_dir/Trinity.fasta
~~~

Which will generate the following data:

~~~ {.output}
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':    377
Total trinity transcripts:    382
Percent GC: 39.25

########################################
Stats based on ALL transcript contigs:
########################################

    Contig N10: 2576
    Contig N20: 1917
    Contig N30: 1594
    Contig N40: 1377
    Contig N50: 1123

    Median contig length: 429
    Average contig: 709.04
    Total assembled bases: 270853


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

    Contig N10: 2532
    Contig N20: 1913
    Contig N30: 1579
    Contig N40: 1365
    Contig N50: 1075

    Median contig length: 427
    Average contig: 697.11
    Total assembled bases: 262812
~~~

This summary tells us:

* How many genes and transcripts were assembled
* GC content
* Statistics on the average size of the contigs
* Number of bases assembled
* Statistics on the average size of the longest isoform of each gene

From this summary, the concept of Contig N50, N40, etc. is particularly important.
For example, N50 indicates the size of the medium contig (or 50%) when
all contigs are sorted by size. N40 is 40% etc. This measure can help to not simply rely on the longest contig and also allows us to observe
how the overall length is increasing in all the transcripts.

Finally, we will perform a blast of the first 5 sequences
to identify which organism they come from. Navigate to:

[NCBI Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi)

> ## Challenge - Assemble with all libraries {.challenge}
>
> Now try assembling all samples together to create a reference de novo transcriptome
>  

> ## What's next? {.callout}
> 
> There are several things that can be done with this de novo reference:
>
> * Check how complete the assembly is using [BUSCO](https://busco.ezlab.org/)
> * Annotate the transcripts
> * Carry out differential gene expression analysis using transcripts as a reference
>


