nadia-reads
===========

Overview
--------

Perform quality control on sequencing reads and trim adapter sequences.

The purpose of this step is checking the quality of sequencing reads before doing
the alignment. The reads from single cell experiments could contain adapter
sequences as well as polyA/polyT sequences if they are derived from short RNA 
molecules. Those sequences need to be trimmed out from the reads to increase the 
mapping rate to reference genome.


How it works
------------

1. Concatenate reads across lanes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequencing reads for one sample could be delivered in muliple fastq files for
different lanes. So first, they will be concatenated into a single fastq file.



.. tip:: 

    We can specify multiple fastq file in `-r1, --read1`_ and `-r2, --read2`_ 
    arguments, but they must have the same order. For example:

    .. code-block:: bash

        nadia-reads \
            -r1 lane1_R1.fastq.gz lane2_R1.fastq.gz \
            -r1 lane1_R2.fastq.gz lane2_R2.fastq.gz \

.. note:: 

    When multiple fastq files are input, sample name is required by `-n, --name`_
    argument.


2. Trim sequences by Cutadapt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To trim sequences, use `--trim`_ flag.

The sequences to trim could be specified by `-a, --adapter`_ argument. This
argument accepts an FASTA file containing trim sequences.

.. _default-adapter-sequence:

If no adapter file is specified, the following sequences will be used by default:

.. code-block::

    >Illumina_Universal
    AGATCGGAAGAG
    >PrefixNX/1
    AGATGTGTATAAGAGACAG
    >Trans1
    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
    >Trans1_rc
    CTGTCTCTTATACACATCTGACGCTGCCGACGA
    >Trans2
    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
    >Trans2_rc
    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
    >polyA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    >polyT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    >polyC
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    >polyG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    >drop-seq
    GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGC
    >Nextera 
    CTGTCTCTTATACACATCT


Besides that, Cutadapt is run with the following default parameters. 
See `cutadapt documentation <https://cutadapt.readthedocs.io/en/stable/>`_ 
for more detail.

.. code-block:: bash

    cutadapt \
        --max-n=0 \
        --minimum-length=20 \
        -q 20,20 \
        --overlap=8


1. Quality control by FASTQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quality of sequencing reads is reported by FASTQC. FASTQC is run on both raw
reads and trimmed reads. Adapter contamination is checked agaist the same 
sequences as Cutadapt (specified by `-a, --adapter`_ argument).


Input
-----

* Sequencing reads in FASTQ format
* Adapter sequences in FASTA format



Output
------

* A multiqc report of FASTQC and Cutadapt. You can download an 
  :download:`example report </_static/reports/nadia_reads_report.html>` 
* Trimmed reads in FASTQ format (ready to be aligned)


Usage examples
--------------

.. code-block:: bash

    nadia-reads \
        -r1 tests/testdata/L1_R1.fastq.gz tests/testdata/L2_R1.fastq.gz \
        -r2 tests/testdata/L1_R2.fastq.gz tests/testdata/L2_R2.fastq.gz \
        -n test_sample \
        -o tests/testresult/reads \
        --trim -a tests/testdata/adapters.fa




Argument details
----------------

``-r1``, ``--read1``
~~~~~~~~~~~~~~~~~~~~
*Required*

Read 1 fastq files. If multiple files are input, they must have the same order with --read2

``-r2``, ``--read2``
~~~~~~~~~~~~~~~~~~~~
*Required*

Read 1 fastq files. If multiple files are input, they must have the same order with --read2

``-n``, ``--name``
~~~~~~~~~~~~~~~~~~~~

Sample name. It will be used for naming output files.

Required if there are multiple input files

If single fastq file is input and **--name** is not specified, then filename of
read 2 will be used for sample name.


``-o``, ``--outdir``
~~~~~~~~~~~~~~~~~~~~
*Required*

Output directory.

``--trim``
~~~~~~~~~~

If this flag is used, then run Cutadapt to trim sequence


``-a``, ``--adapter``
~~~~~~~~~~~~~~~~~~~~~

Path to adapter file in FASTA format. See `2. Trim sequences by Cutadapt`_
