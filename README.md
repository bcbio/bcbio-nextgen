![bcbio banner](https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/docs/contents/images/banner.png)

Validated, scalable, community developed variant calling, RNA-seq and small RNA analysis. You write a high level configuration file specifying your inputs and analysis parameters. This input drives a parallel run that handles distributed execution, idempotent processing restarts and safe transactional steps. bcbio provides a shared community resource that handles the data processing component of sequencing analysis, providing researchers with more time to focus on the downstream biology.

[![Build Status](https://travis-ci.org/bcbio/bcbio-nextgen.svg?branch=master)](https://travis-ci.org/bcbio/bcbio-nextgen)
[![Documentation status](https://readthedocs.org/projects/bcbio-nextgen/badge/?version=latest)](https://bcbio-nextgen.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3564938.svg)](https://doi.org/10.5281/zenodo.3564938)

## Features

* Community developed: We welcome contributors with the goal of overcoming the biological, algorithmic and computational challenges that face individual developers working on complex pipelines in quickly changing research areas. See our [users page](https://bcbio-nextgen.readthedocs.io/en/latest/contents/about.html#users) for examples of bcbio-nextgen deployments, and the [developer documentation](https://bcbio-nextgen.readthedocs.io/en/latest/contents/development.html) for tips on contributing.
* Installation: [A single installer script](https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#automated) prepares all third party software, data libraries and system configuration files.
* [Automated validation](https://bcb.io/2014/05/12/wgs-trio-variant-evaluation/): Compare variant calls against common reference materials or sample specific SNP arrays to ensure call correctness. Incorporation of multiple approaches for alignment, preparation and variant calling enable unbiased comparisons of algorithms.
* Distributed: Focus on [parallel analysis and scaling](https://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/) to handle large population studies and whole genome analysis. Runs on single multicore computers, in compute clusters using [IPython parallel](https://ipyparallel.readthedocs.io/en/latest/), or on the Amazon cloud. See the [parallel documentation](https://bcbio-nextgen.readthedocs.org/en/latest/contents/parallel.html) for full details.
* Multiple analysis algorithms: bcbio-nextgen provides configurable [variant calling, RNA-seq and small RNA pipelines](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html).

## Quick start

1. [Install](https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#automated) `bcbio-nextgen` with all tool dependencies and data files:
    ```shell script
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py /usr/local/share/bcbio --tooldir=/usr/local \
          --genomes GRCh37 --aligners bwa --aligners bowtie2
    ```
   producing an editable [system configuration file](https://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_system.yaml) referencing the installed software, data and system information.

1. [Automatically create a processing description](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#automated-sample-configuration) of sample FASTQ and BAM files from your project, and a CSV file of sample metadata:
    ```shell script
    bcbio_nextgen.py -w template freebayes-variant project1.csv sample1.bam sample2_1.fq sample2_2.fq
    ```
    This produces a [sample description file](https://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_sample.yaml) containing pipeline [configuration options](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html).

1.  Run analysis, distributed across 8 local cores:
    ```shell script
    cd project1/work
    bcbio_nextgen.py ../config/project1.yaml -n 8
    ```

## Documentation

See the [full documentation](https://bcbio-nextgen.readthedocs.io/en/latest/) and [longer analysis-based articles](https://bcb.io). We welcome enhancements or problem reports using [GitHub](https://github.com/bcbio/bcbio-nextgen/issues) and discussion on the [biovalidation mailing list](https://groups.google.com/d/forum/biovalidation).

## Contributors

* [Miika Ahdesmaki](https://github.com/mjafin), AstraZeneca
* [Luca Beltrame](https://github.com/lbeltrame), IRCCS "Mario Negri" Institute for Pharmacological Research, Milan, Italy
* [Christian Brueffer](https://github.com/cbrueffer), Lund University, Lund, Sweden
* [Alla Bushoy](https://github.com/abushoy), AstraZeneca
* [Guillermo Carrasco](https://github.com/guillermo-carrasco), Science for Life Laboratory, Stockholm
* [Nick Carriero](https://www.simonsfoundation.org/team/nick-carriero/), Simons Foundation
* [Brad Chapman](https://github.com/chapmanb), Harvard Chan Bioinformatics Core
* [Saket Choudhary](https://github.com/saketkc), University Of Southern California
* [Peter Cock](https://github.com/peterjc), The James Hutton Institute
* [Matthias De Smet](https://github.com/matthdsm), Center for Medical Genetics, Ghent University Hospital, Belgium
* [Matt Edwards](https://github.com/matted), MIT
* [Mario Giovacchini](https://github.com/mariogiov), Science for Life Laboratory, Stockholm
* [Karl Gutwin](https://twitter.com/kgutwin), Biogen
* [Jeff Hammerbacher](https://github.com/hammer), Icahn School of Medicine at Mount Sinai
* [Oliver Hofmann](https://umccr.github.io/), University of Melbourne Centre for Cancer Research
* [John Kern](https://github.com/kern3020)
* [Rory Kirchner](https://github.com/roryk), Harvard Chan Bioinformatics Core
* [Tetiana Khotiainsteva](https://github.com/tetianakh), Ardigen
* [Kerrin Mendler](https://github.com/kmendler), AstraZeneca
* [Sergey Naumenko](https://github.com/naumenko-sa), Harvard Chan Bioinformatics Core
* [Jakub Nowacki](https://github.com/jsnowacki), AstraZeneca
* [John Morrissey](https://github.com/jwm), Harvard Chan Bioinformatics Core
* [Lorena Pantano](https://github.com/lpantano), Harvard Chan Bioinformatics Core
* [Brent Pedersen](https://github.com/brentp), University of Colorado Denver
* [James Porter](https://github.com/porterjamesj), The University of Chicago
* [Vlad Saveliev](https://github.com/vladsaveliev), Center for Algorithmic Biotechnology, St. Petersburg University
* [Valentine Svensson](https://github.com/vals), Science for Life Laboratory, Stockholm
* [Paul Tang](https://github.com/tanglingfung), UCSF
* [Stephen Turner](https://github.com/stephenturner), University of Virginia
* [Roman Valls](https://github.com/brainstorm), Science for Life Laboratory, Stockholm
* [Kevin Ying](https://github.com/kevyin), Garvan Institute of Medical Research, Sydney, Australia

## License

The code is freely available under the [MIT license](https://opensource.org/licenses/mit-license.html).
