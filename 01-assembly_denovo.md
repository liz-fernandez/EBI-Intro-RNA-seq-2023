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
$ cd RNAseq-2023/
$ mkdir De_Novo_Assembly
$ cd De_Novo_Assembly
~~~ 

And we will copy the fastq files to our directory:

~~~ {.bash}
$ cp /home/training/RNAseq-2023/DATA/RNASEQ_data.zip .
~~~ 

And decompress it:

~~~ {.bash}
$ unzip RNASEQ_data.zip
$ cd RNASEQ_data
~~~ 

And check everything is as expected:

~~~ {.bash}
$ ls Sp*fq.gz 
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

And use Trinity with a generic command:

~~~ {.bash}
$ Trinity --seqType fq \
--left Sp_log.left.fq.gz,Sp_hs.left.fq.gz \
--right Sp_log.right.fq.gz,Sp_hs.right.fq.gz \
--CPU 2 --max_memory 2G --min_contig_length 150
~~~

This command will take approximately 15 minutes to assemble a
transcriptome.

The options (flags) that we have used are the following:

* --seqType fq - We indicate that we are using fastq type files
* --left - Left side reads (o R1)
* --right - Right side reads (o R2)
* --CPU 2 - Use 2 CPUs
* --max_memory 2G - Indicate that the maximum RAM to use is 2GB
* --min_contig_lenght - Minimum size for contigs (default is 200)

These are the most essential options to carry out the analysis.

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

     ______  ____   ____  ____   ____  ______  __ __
    |      ||    \ |    ||    \ |    ||      ||  |  |
    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
      |__|  |__|\_||____||__|__||____|  |__|  |____/

    Trinity-v2.15.0


#
#
# Required:
#
#  --seqType <string>      :type of reads: ('fa' or 'fq')
#
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provided in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
#
#  Or,
#      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.
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
#                                   (def=200, must be >= 100)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#                                   (** note: experimental parameter **, this functionality continues to be under development)
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --long_reads_bam <string>       :long reads to include for genome-guided Trinity
#                                  (bam file consists of error-corrected or circular consensus (CCS) pac bio read aligned to the genome)
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
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/home/training/trinity_out_dir" 
#                                    note: must include 'trinity' in the name as a safety precaution! )
#  
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (Trinity-v2.15.0) and exits.
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
#            (if you have multiple samples, use --samples_file ... see above for details)
#
#    and for Genome-guided Trinity, provide a coordinate-sorted bam:
#
#        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
#
#     see: /usr/local/bin/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#
#     For more details, visit: http://trinityrnaseq.github.io
#
###############################################################################
~~~ 

Our work must be finished, let's review the transcripts
generated, which are in the file `Trinity.fasta`:

~~~ {.bash}
$ head trinity_out_dir.Trinity.fasta
~~~

We note that the results are sequences in `fasta` format.

~~~ {.output}
>TRINITY_DN83_c0_g1_i1 len=400 path=[0:0-399]
TGTTGGGAATAATTCTGTCATACCATACCGTGAATGTCGCAGAATGTTTCAATTTCAGTA
GCACTGGCATTTTATCTTTTCTGACTTATTTTCAGTCCGGCTGTTCCCTTCTCTTTTGGT
GTCCCTGTTGTCTTGGACAATTTGTCACCCCAAGGACCCTTGTACTTGGTTGTTTAAACA
ATCTCCGTAACCTCATTTCTTTTCTCCTGTCTCTTGCATTTCACACACCAATCCAAACAA
GTTCAGCGTGAGGACTGCTTTAACATCCGTGATCGAGCCTAGACGGTTGTCCCAATGTTT
AACAAGGAGCTCGATGTACTTTTTTATTTGGTCCGTTGACTTTCCGCAAACATGCAGTGA
GAGAATACGGAAGTGAAGAGAATTTCGGTAATGTATCGCA
>TRINITY_DN50_c0_g1_i1 len=182 path=[0:0-181]
ATCATATGGGTCATTTACGAACCCGAGATGGATCCACTCATTCAACATTTTTGTATAAGG
~~~

## Analyzing the statistics of the assembled transcriptome

Let's count the number of assembled transcripts:

~~~ {.bash}
$ grep '>' trinity_out_dir.Trinity.fasta | wc -l
~~~

~~~ {.output}
439
~~~

How many transcripts were assembled? This is somewhat informative, but be aware
that, the deeper you sequence, the more transcripts you will reconstruct. When
doing a *de novo* assembly it is common to assemble substantially more transcripts
than what you would expect when doing a transcriptome assembly guided by a
reference genome. 

We can capture some additional statistics about this assembly
using a program that is part of Trinity:

~~~ {.bash}
$ TrinityStats.pl trinity_out_dir.Trinity.fasta
~~~

Which will generate the following data:

~~~ {.output}
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	431
Total trinity transcripts:	439
Percent GC: 38.97

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 2576
	Contig N20: 2294
	Contig N30: 1774
	Contig N40: 1423
	Contig N50: 1152

	Median contig length: 398
	Average contig: 678.84
	Total assembled bases: 298012


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 2576
	Contig N20: 2190
	Contig N30: 1670
	Contig N40: 1388
	Contig N50: 1143

	Median contig length: 385
	Average contig: 663.88
	Total assembled bases: 286133

~~~

This summary tells us:

* How many genes and transcripts were assembled (should be the same number as we
  calculated before using 'grep')
* GC content
* Statistics on the average size of the contigs
* Number of bases assembled
* Statistics on the average size of the longest isoform of each gene

From this summary, the concept of Contig N50, N40, etc. is particularly important.
For example, N50 indicates that at least half of the assembled bases are in
contigs of at least that contig lenght. N40 is 40% etc. This measure can help to not simply rely on the longest contig and also allows us to observe
how the overall length is increasing in all the transcripts.

However, summary statistics of the assembly are not very informative. 

A useful metric is to assess the number of fully reconstructed coding transcripts.
This can be calculated by aligning all transcripts to a database of coding
sequences (such as [UNIPROT](uniprot.org)) with an aligner (such as BLASTX). We are going to a mock search against a smaller database for time's sake:

Downloading and decompressiong short database:

~~~ {.bash}
$ wget https://github.com/liz-fernandez/EBI-Intro-RNA-seq-2023/raw/main/DATA/mini_sprot.pep.gz
$ gunzip mini_sprot.pep.gz
~~~ 

Building blast database:

~~~ {.bash}
makeblastdb -in mini_sprot.pep -dbtype prot -parse_seqids -out mini_sprot_pep.db
~~~

~~~ {.output}
Building a new DB, current time: 02/21/2023 12:33:57
New DB name:   /home/training/RNAseq-2023/De_Novo_Assembly/RNASEQ_data/mini_sprot_pep.db
New DB title:  mini_sprot.pep
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 408 sequences in 0.008569 seconds.
~~~

Doing blast search against generated database:

~~~ {.bash}
$ blastx -query trinity_out_dir.Trinity.fasta \
	-db mini_sprot_pep.db -out blastx.outfmt6 \
	-evalue 1e-20 -num_threads 2 -max_target_seqs 1 -outfmt 6
~~~

Then we use a Trinity script to calculate the overall size:

~~~ {.bash}
$ analyze_blastPlus_topHit_coverage.pl \
	blastx.outfmt6 trinity_out_dir.Trinity.fasta \
	mini_sprot.pep | column -t
~~~

~~~ {.output}
#hit_pct_cov_bin  count_in_bin  >bin_below
100               10            10
90                3             13
80                2             15
70                3             18
60                1             19
50                1             20
40                0             20
30                0             20
20                1             21
10                0             21
~~~

The table above lists the number of proteins that are matched by percentage of
sequenced map in bins (2nd column). For example, 10 proteins have between 90% and 100% of their lenght aligned to the database, 3 between 80% and 90% and so on. The third column is a running total, just summing all transcripts that are above that percentage threshold. 

This type of mapping provides and estimate to whether you have sequenced deep
enough or you need to do more RNA-seq to assemble a more complete transcriptome. 

A more low througput way to do check if the assembled sequences make any sense is
to a blast of a few sequences, to identify which organism they come from. Navigate to:

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


