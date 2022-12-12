System requirements
===================

System
------
Nadiatools requires a linux system. It is extensively tested on Ubuntu.

STARsolo pipeline requires more RAM and computational power than Alevin-fry 
pipeline. It is recommended to have at least 8-core CPUs and 32GB of RAM. The
disk space depends on the amount of data (500GB recommended).

Dependencies
------------
Nadiatools require following softwares to be installed:

.. csv-table::
    :header: Softwares,Version
    :widths: 100,50

    fastqc,0.11.9
    cutadapt,4.1
    star,2.7.10a
    salmon,1.9.0
    alevin-fry,0.8.0
    AlevinQC,1.14.0
    multiqc,1.13
    anndata,0.8.0
    scanpy,1.9.1
    scipy,1.9.3
    harmonypy,0.0.6
    scanorama,1.7.3
    scrublet,0.2.3
    umap-learn,0.5.3
    matplotlib,3.6.1
    numpy,1.23.4
    pandas,1.5.1
    plotly,5.11.0
    seaborn,0.11.2
    biopython,1.79
    sphinx,5.3.0
    sphinx-rtd-theme,1.1.1
    python,3.7
    R,4.2.2