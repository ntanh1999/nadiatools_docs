Overview
========

Nadiatools performs *pre-processing* and *data processing* steps in 
`a scRNA-seq analysis workflow <http://bioconductor.org/books/3.16/OSCA.intro/analysis-overview.html#outline>`_. 
It is specially designed to analyze data from Nadia instrument, which is based 
on Drop-seq technology.

What we can do
--------------

* Quality control on sequencing reads
* Trim adapter sequences and low quality bases
* Perform read alignment and generate gene expression matrix by STARsolo or Alevin-fry
* Support single-cell and single-nucleus workflow
* Perform data processing steps: (1) Calculating quality controls metrics, (2)
  Doublet detection, (3) Normalizing data, (4) Feature selection, (5) Dimentionality 
  reduction, (6) Visualization (UMAP, T-SNE), (7) Integrating datasets.

Noticeable features
-------------------

* Designed for Nadia instrument:
    * Barcode and UMI structures are built-in
    * No barcode whitelist is required
* Produce an intensive report in html format
* Produce many useful plots (knee plot, violin plot, PCA, UMAP, T-SNE plots,...)
* Output files could be used for further downstream analyses (e.g. AnnData 
  object for cellxgene visualization tool, mtx matrix for ASAP tool,...)


Command description
-------------------

* :doc:`/Commands/ref`: prepare index for alignment. It can generate both STAR 
  index and `splici-index <https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/>`_ 
  for Salmon-Alevinfry workflow.
* :doc:`/Commands/reads`: concatenate reads files across lanes, perform read 
  quality control, trim adapter sequences and low quality bases.
* :doc:`/Commands/quant`: takes reads files from ``nadia-reads`` and 
  index from ``nadia-ref`` and performs alignment, gene 
  expression quantification and barcode filtering. It supports both STARsolo 
  and Alevin-fry pipeline, single-cell and single-nucleus workflow. 
  ``nadia-quant`` produce a feature-barcode matrix and some plots for 
  further processing.
* :doc:`/Commands/process`: takes a matrix (mtx or h5ad) from ``nadia-quant`` 
  and performs *data processing* steps: calculate qc metrics, filter cells, 
  filter genes, doublet detection, normalizing, feature selection, 
  dimentionality reduction and visualization.
* :doc:`/Commands/combine`: takes processed h5ad files from ``nadia-process``
  and integrates multiple sample together.