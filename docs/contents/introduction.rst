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

- `Science for Life Laboratory`_: The genomics core platform in
  the Swedish National Infrastructure (`NGI`_) for genomics, has crunched
  over 16TBp (terabasepairs) and processed almost 7000+ samples
  from the beggining of 2013 until the end of July. `UPPMAX`_, our
  cluster located in Uppsala runs the pipeline in production since 2010.

.. _Massachusetts General Hospital: http://molbio.mgh.harvard.edu/
.. _Galaxy: http://galaxyproject.org/

- `Science for Life Laboratory`_: The genomics core platform in
  the Swedish National Infrastructure (`NGI`_) for genomics, has crunched
  over 16TBp (terabasepairs) and processed almost 7000+ samples
  from the beggining of 2013 until the end of July. `UPPMAX`_, our
  cluster located in Uppsala runs the pipeline in production since 2010.

.. _Science for Life Laboratory: http://scilifelab.se/
.. _NGI: https://portal.scilifelab.se/genomics/
.. _UPPMAX: http://www.uppmax.uu.se/uppnex

- `Institute of Human Genetics, UCSF`_: The Genomics Core Facility
  utilizes bcbio-nextgen in processing more than 2000 whole genome,
  exome, RNA-seq, ChIP-seq on various projects. This pipeline
  tremendously lowers the barrier of getting access to next generation
  sequencing technology. The community engaged here is also very
  helpful in providing best practices advices and up-to-date solution
  to ease scientific discovery.

.. _Institute of Human Genetics, UCSF: http://humangenetics.ucsf.edu/

- `IRCCS "Mario Negri" Institute for Pharmacological Research`_:
  The Translational Genomics Unit in the Department of Oncology uses
  bcbio-nextgen for targeted resequencing (using an Illumina MiSeq) to
  identify mutations and other variants in tumor samples to
  investigate their link to tumor progression, patient survival and
  drug sensitivity and resistance. A
  `poster from the 2014 European Society of Human Genetics meeting`_
  provides more details on usage in ovarian cancer.

.. _IRCCS "Mario Negri" Institute for Pharmacological Research: http://www.marionegri.it
.. _poster from the 2014 European Society of Human Genetics meeting: https://github.com/chapmanb/bcbb/raw/master/posters/beltrame_ESHG_poster_05_2014.reduced.pdf

- `The Translational Genomics Research Institute (TGen)`_: 
  Members of the `Huentelman lab`_ at TGen apply bcbio-nextgen to a wide 
  variety of studies of with a major focus in the neurobiology of aging 
  and neurodegeneration in collaboration with the The Arizona Alzheimer's Consortium (`AAC`_)
  and  the `McKnight Brain Research Foundation`_.
  We also use bcbio in studies of rare diseases in children through TGen's 
  Center for Rare Childhood Disorders (`C4RCD`_),  and other rare diseases such as 
  Multiple System Atrophy (`MSA`_). bcbio-nextgen has also been instrumental in 
  projects for TGen's Program for Canine Health & Performance (`PCHP`_) 
  and numerous RNA-seq projects using rodent models. Our work with bcbio 
  started with a parnership with `Dell` and The Neuroblastoma and 
  Medulloblastoma Translational Research Consortium (`NMTRC`_), 
  and TGen as part of a Phase I clinical trial in these rare childhood cancers.

.. _The Translational Genomics Research Institute (TGen): http://www.tgen.org
.. _Huentelman lab: http://www.tgen.org/research/research-faculty/matt-huentelman.aspx
.. _AAC: http://www.azalz.org
.. _McKnight Brain Research Foundation: http://tmbrf.org  
.. _C4RCD: http://www.c4rcd.org
.. _MSA: http://www.tgen.org/research/multiple-system-atrophy-(msa)-research-registry.aspx
.. _PCHP: http://www.tgen.org/research/canine-health-performance.aspx
.. _Dell: http://www.dell.com/learn/us/en/70/healthcare
.. _NMTRC: http://nmtrc.org/about
   
