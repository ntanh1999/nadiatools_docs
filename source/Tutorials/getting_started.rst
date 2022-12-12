Getting started
=========================

Introduction
------------
This tutorial go through steps of analyzing a scRNA-seq sample by nadiatools, 
from create index, reads quality control, alignment to visualization. 

About the example data:

* Sample name: HEK_3T3
* Use RNAdia2.0 reagent kit
* mixed human 293T and mouse 3T3 cells
* read 1 is 28 bp, read 2 is 91 bp

Running this tutorial using Alevin-fry pipeline requires at least 16 GB of memory. 
With STARsolo pipeline, the memory requirement is much larger (57 GB if not create STAR index, 122 GB otherwise)

Download sample data
--------------------

First, we download fastq files into ``RAWDATA`` folder. 
You can download `here <https://figshare.com/articles/dataset/HEK_3T3_samle/21709907>`_.

.. code-block:: bash

    mkdir RAWDATA
    cd RAWDATA
    wget https://figshare.com/ndownloader/files/38515226 -O 1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz
    wget https://figshare.com/ndownloader/files/38515232 -O 1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz
    wget https://figshare.com/ndownloader/files/38515238 -O 2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz
    wget https://figshare.com/ndownloader/files/38515244 -O 2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz



.. code-block:: bash

    .
    └── RAWDATA
        ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz
        ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz
        ├── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz
        └── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz


Prepare reference
-----------------
This mixed species experiment requires a hybrid genome of human and mouse. 
Reference genome and transcript annoation could be downloaded from this 
`link <https://figshare.com/account/collections/6340910>`_.
We will download ``homo_38.106_mus_39.106.fa.gz`` and ``homo_38.106_mus_39.106.gtf.gz`` in to ``REF`` folder.

.. code-block:: bash

    mkdir REF
    cd REF
    wget https://figshare.com/ndownloader/files/38514659 -O homo_38.106_mus_39.106.gtf.gz
    wget https://figshare.com/ndownloader/files/38514680 -O homo_38.106_mus_39.106.fa.gz



.. code-block:: bash

    .
    │── RAWDATA
    │    ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz
    │    ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz
    │    ├── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz
    │    └── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz
    └── REF
        ├── homo_38.106_mus_39.106.fa.gz
        └── homo_38.106_mus_39.106.gtf.gz


We are now using ``nadia-ref`` to generate index for read alignment. In this tutorial, 
we use STAR as aligner, so we will create a STAR index. We also perform filtering 
on GTF file and only keep "protein_coding" gene biotype. 

.. code-block:: bash

    nadia-ref \
        -g REF/homo_38.106_mus_39.106.fa.gz \
        -a REF/homo_38.106_mus_39.106.gtf.gz \
        -i star \
        -l 91 \
        -o REF \
        -f  \
        --gene-biotype protein_coding


By default, we used all CPUs (``-t 0``) to create index. In our case, all 24 cores are used. 
It took us 34 minutes and required 122 GB of memory to create this STAR index. 
The STAR index folder itself requires 56.3 GB of disk storage. If your computer
cannot handle STAR index, you can download a pre-made one (https://ewels.github.io/AWS-iGenomes/) 
or use Salmon as an aligner.

.. Command being timed: "nadia-ref -g REF/homo_38.106_mus_39.106.fa.gz -a REF/homo_38.106_mus_39.106.gtf.gz -i star -l 91 -o REF -f --gene-biotype protein_coding"
.. User time (seconds): 25589.14
.. System time (seconds): 341.14
.. Percent of CPU this job got: 1269%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 34:03.06
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 128256572
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 3564233
.. Minor (reclaiming a frame) page faults: 106855140
.. Voluntary context switches: 3584421
.. Involuntary context switches: 3700039
.. Swaps: 0
.. File system inputs: 156091288
.. File system outputs: 305201016
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0


Use the following command to create salmon index. It took us 18 minutes and 13.8 
GB of memory. The size of salmon index folder is 16.5 GB.

.. code-block:: bash

    nadia-ref \
        -g REF/homo_38.106_mus_39.106.fa.gz \
        -a REF/homo_38.106_mus_39.106.gtf.gz \
        -i salmon \
        -l 91 \
        -o REF \
        -f  \
        --gene-biotype protein_coding    

.. Command being timed: "nadia-ref -g REF/homo_38.106_mus_39.106.fa.gz -a REF/homo_38.106_mus_39.106.gtf.gz -i salmon -l 91 -o REF -f --gene-biotype protein_coding"
.. User time (seconds): 10979.14
.. System time (seconds): 69.67
.. Percent of CPU this job got: 1003%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 18:21.42
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 14493240
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 258214
.. Minor (reclaiming a frame) page faults: 25739781
.. Voluntary context switches: 388029
.. Involuntary context switches: 1469988
.. Swaps: 0
.. File system inputs: 3543344
.. File system outputs: 84253080
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0


Here is the output:

.. code-block:: bash

    .
    ├── RAWDATA
    │   ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz
    │   ├── 1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz
    │   ├── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz
    │   └── 2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz
    └── REF
        ├── homo_38.106_mus_39.106.fa.gz
        ├── homo_38.106_mus_39.106.gtf.gz
        ├── homo_38.106_mus_39.106_filtered.gtf
        ├── salmon_index
        └── star_index



Reads quality control
---------------------

In this step, we perform quality control on sequencing reads and trim adapter sequences by ``nadia-reads``.

Sequencing reads for this sample are delivered in multiple fastq files for 
different lanes. We can specify multiple fastq file in ``-r1, --read1`` and ``-r2, --read2``
arguments, but they must have the same order.

For this sample, we will trim adapter sequences by cutadapt and use :ref:`the default
adapter sequences <default-adapter-sequence>`.

.. code-block:: bash

    nadia-reads \
        -r1 RAWDATA/1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz RAWDATA/2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz \
        -r2 RAWDATA/1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz RAWDATA/2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz \
        -n HEK_3T3 \
        -o results \
        --trim


.. Command being timed: "nadia-reads -r1 RAWDATA/1_HFKGJDMXY.1.ACCTTCTC.unmapped.1.fastq.gz RAWDATA/2_HFKGJDMXY.2.ACCTTCTC.unmapped.1.fastq.gz -r2 RAWDATA/1_HFKGJDMXY.1.ACCTTCTC.unmapped.2.fastq.gz RAWDATA/2_HFKGJDMXY.2.ACCTTCTC.unmapped.2.fastq.gz -n HEK_3T3 -o results/ReadQC --trim"
.. User time (seconds): 10526.16
.. System time (seconds): 170.81
.. Percent of CPU this job got: 682%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 26:08.23
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 4090208
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 1344
.. Minor (reclaiming a frame) page faults: 37241982
.. Voluntary context switches: 3689886
.. Involuntary context switches: 3964794
.. Swaps: 0
.. File system inputs: 20598296
.. File system outputs: 38503328
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0


Here is the output files:

.. code-block:: bash

    .
    ├── adapter.tsv
    ├── fastqc
    ├── nadia_reads_report.html
    ├── raw
    │   ├── HEK_3T3_R1.fastq.gz
    │   └── HEK_3T3_R2.fastq.gz
    └── trimmed
        ├── HEK_3T3_cutadapt.log
        ├── HEK_3T3_trimmed_R1.fastq.gz
        └── HEK_3T3_trimmed_R2.fastq.gz


You can download :download:`nadia_reads_report.html </_static/reports/nadia_reads_report.html>`

``HEK_3T3_trimmed_R1.fastq.gz`` and ``HEK_3T3_trimmed_R1.fastq.gz`` are the 
trimmed reads, which will be used in the next step.

Align and quanfify gene expression
----------------------------------

Giving the trimmed reads and the contructed index, the next step is align reads
to genome and quantify gene expression. 

The next command is to run STARsolo pipeline:

.. code-block:: bash

    nadia-quant \
        -r1 results/trimmed/HEK_3T3_trimmed_R1.fastq.gz \
        -r2 results/trimmed/HEK_3T3_trimmed_R2.fastq.gz \
        -i REF/star_index \
        -o results \
        -w single-cell \
        -a starsolo \
        -s RNAdia \
        --top-cells 2500 \
        --mito "^MT-" --ribo "^RP[SL]"


Explaination:

* ``-r1 results/trimmed/HEK_3T3_trimmed_R1.fastq.gz``, ``-r2 results/trimmed/HEK_3T3_trimmed_R2.fastq.gz``: the trimmed reads from ``nadia-reads``
* ``-i REF/star_index``: star index folder from ``nadia-ref``
* ``-o results``: output to results folder.
* ``-w single-cell``: use single-cell analysis workflow. (for single-nucleus sample, use ``-w single-nucleus`` )
* ``-a starsolo``: use STARsolo pipeline
* ``-s RNAdia``: this sample use RNAdia 2.0 reagent kit, so we will use :ref:` RNAdia barcode structure <rnadia-barcode-structure>`.
* ``--top-cells 2500``: keep 2500 top barcodes
* ``--mito "^MT-" --ribo "^RP[SL]"``: these arguments provides the regular expression for mitochondrial genes and ribosomal genes. 
  They are used to calculate the percentage of these genes per cell.

It took us 10 minutes and 57 GB of memory (using 24 CPUs). 

.. Command being timed: "nadia-quant -r1 results/ReadQC/trimmed/HEK_3T3_trimmed_R1.fastq.gz -r2 results/ReadQC/trimmed/HEK_3T3_trimmed_R2.fastq.gz -i REF/star_index -o results -w single-cell -a starsolo -s RNAdia --top-cells 2500 --mito ^MT- --ribo ^RP[SL]"
.. User time (seconds): 7558.06
.. System time (seconds): 156.98
.. Percent of CPU this job got: 1278%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 10:03.44
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 59971324
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 499
.. Minor (reclaiming a frame) page faults: 32251611
.. Voluntary context switches: 868256
.. Involuntary context switches: 1797211
.. Swaps: 0
.. File system inputs: 50170920
.. File system outputs: 193330720
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0

Here is the output:

.. code-block:: bash

   results
    ├── anndata
    │   ├── HEK_3T3_filter.h5ad
    │   └── HEK_3T3_raw.h5ad
    ├── MTX
    │   └── HEK_3T3
    │       ├── filter
    │       │   ├── barcodes.tsv.gz
    │       │   ├── features.tsv.gz
    │       │   └── matrix.mtx.gz
    │       └── raw
    │           ├── barcodes.tsv.gz
    │           ├── features.tsv.gz
    │           └── matrix.mtx.gz
    ├── nadia_quant_report.html
    └── STARsolo
        └── HEK_3T3

You can download :download:`nadia_quant_report.html </_static/reports/nadia_quant_report.html>`

With STARsolo pipeline, two matrix will be output (``raw`` and ``filter``). Raw 
matrix is created without cell filtering step. On the other hand, the filtered matrix, 
in this case, only contains top 2500 barcodes.

``STARsolo`` folder contains other output files of STARsolo.

To run Alevin-fry pipeline, run the following command. It took 4 minutes and 14 GB of memory.

.. code-block:: bash

    nadia-quant \
        -r1 results/trimmed/HEK_3T3_trimmed_R1.fastq.gz \
        -r2 results/trimmed/HEK_3T3_trimmed_R2.fastq.gz \
        -i REF/salmon_index \
        -o results_2 \
        -w single-cell \
        -a alevin-fry \
        -s RNAdia \
        --top-cells 2500 \
        --mito "^MT-" --ribo "^RP[SL]"


You can download :download:`nadia_quant_alevinfry_report.html </_static/reports/nadia_quant_alevinfry_report.html>`

.. User time (seconds): 1153.99
.. System time (seconds): 20.07
.. Percent of CPU this job got: 524%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 3:43.66
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 14980420
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 2539536
.. Minor (reclaiming a frame) page faults: 2140894
.. Voluntary context switches: 5309053
.. Involuntary context switches: 101113
.. Swaps: 0
.. File system inputs: 27262560
.. File system outputs: 9658704
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0

Process gene expression matrix
------------------------------

In the final step, we will process the gene expression matrix and create some 
visualization. This can be done with the below command:

.. code-block:: bash

    nadia-process \
        --h5ad results/anndata/HEK_3T3_filter.h5ad \
        -o results \
        --filter-doublet \
        --n-pcs 10 \
        --plot-genes Hsp90b1 Lox


This command takes a filtered matrix in h5ad file. 
It performs cell filtering and gene filtering with all default parameters. 
We also filter doublets with the threshold being automatically selected. We 
only use the first 10 principal components to perform UMAP and T-SNE. Lastly, 
UMAP plot of expression of some genes (Hsp90b1, Lox) are generated.

.. Command being timed: "nadia-process --h5ad results/anndata/HEK_3T3_filter.h5ad -o results --filter-doublet --n-pcs 10 --plot-genes Hsp90b1 Lox LOLC1 SPEN"
.. User time (seconds): 118.71
.. System time (seconds): 16.78
.. Percent of CPU this job got: 408%
.. Elapsed (wall clock) time (h:mm:ss or m:ss): 0:33.18
.. Average shared text size (kbytes): 0
.. Average unshared data size (kbytes): 0
.. Average stack size (kbytes): 0
.. Average total size (kbytes): 0
.. Maximum resident set size (kbytes): 1086996
.. Average resident set size (kbytes): 0
.. Major (requiring I/O) page faults: 343
.. Minor (reclaiming a frame) page faults: 293412
.. Voluntary context switches: 3782
.. Involuntary context switches: 321037
.. Swaps: 0
.. File system inputs: 83680
.. File system outputs: 86704
.. Socket messages sent: 0
.. Socket messages received: 0
.. Signals delivered: 0
.. Page size (bytes): 4096
.. Exit status: 0


Here is the output:

.. code-block:: bash

    .
    ├── nadia_process_report.html
    ├── Plots
    │   └── HEK_3T3
    │       ├── doublet_histogram.png
    │       ├── filter_genes_dispersion.png
    │       ├── highest_expr_genes_filtered.png
    │       ├── pca_loadings.png
    │       ├── pca.png
    │       ├── pca_variance_ratio.png
    │       ├── tsne_metadata.png
    │       ├── umap_doublet.png
    │       ├── umap_genes_of_interest.png
    │       ├── umap_metadata.png
    │       ├── violin_doublet.png
    │       └── violin_metadata_filtered.png
    └── processed_anndata
       └── HEK_3T3_processed.h5ad


``HEK_3T3_processed.h5ad`` could be used directly with cellxgene.

You can download :download:`nadia_process_report.html </_static/reports/nadia_process_report.html>`