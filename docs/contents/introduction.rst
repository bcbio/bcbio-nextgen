Introduction
------------

bcbio-nextgen provides best-practice pipelines for automated analysis
of high throughput sequencing data with the goal of being:

- Quantifiable: Doing good science requires being able to accurately
  assess the quality of results and re-verify approaches as new
  algorithms and software become available.

- Analyzable: Results feed into tools to make it easy to query and
  visualize the results.

- Scalable: Handle large datasets and sample populations on distributed
  heterogeneous compute environments.

- Reproducible: Track configuration, versions, provenance and command
  lines to enable debugging, extensive and reproducibility of results.

- Community developed: The development process is fully open and
  sustained by contributors from multiple institutions. By working
  together on a shared framework, we can overcome the challenges
  associated with maintaining complex pipelines in a rapidly changing
  area of research.

- Accessible: Bioinformaticians, biologists and the general public
  should be able to run these tools on inputs ranging from research
  materials to clinical samples to personal genomes.

Users
=====

A sample of institutions using bcbio-nextgen for solving biological
problems. Please submit your story if you're using the pipeline in
your own research.

- `Harvard School of Public Health`_: We use bcbio-nextgen within the
  bioinformatics core for variant calling on large population studies
  related to human health like Breast Cancer and Alzheimer's disease.
  Increasing scalability of the pipeline has been essential for
  handling study sizes of more than 1400 whole genomes.

.. _Harvard School of Public Health: http://compbio.sph.harvard.edu/chb/

- `Massachusetts General Hospital`_: The Department of Molecular
  Biology uses the pipeline to automatically process samples coming
  off Illumina HiSeq instruments. Automated pipelines perform
  alignment and sample-specific analysis, with results directly
  uploaded into a local `Galaxy`_ instance.

.. _Massachusetts General Hospital: http://molbio.mgh.harvard.edu/
.. _Galaxy: http://galaxyproject.org/

- `Institute of Human Genetics, UCSF`_: The Genomics Core Facility
  utilizes bcbio-nextgen in processing more than 2000 whole genome,
  exome, RNA-seq, ChIP-seq on various projects. This pipeline
  tremendously lowers the barrier of getting access to next generation
  sequencing technology. The community engaged here is also very
  helpful in providing best practices advices and up-to-date solution
  to ease scientific discovery.

.. _Institute of Human Genetics, UCSF: http://humangenetics.ucsf.edu/

- `IRCCS "Mario Negri" Institute for Pharmacological Research`_: The
   Translational Genomics Unit in the Department of Oncology uses
   bcbio-nextgen for targeted resequencing (using an Illumina MiSeq) to
   identify mutations and other variants in tumor samples to investigate their
   link to tumor progression, patient survival and drug
   sensitivity and resistance.

.. _IRCCS "Mario Negri" Institute for Pharmacological Research: http://www.marionegri.it
