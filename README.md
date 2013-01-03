## bcbio-nextgen

A python toolkit providing best-practice pipelines for fully automated high
throughput sequencing analysis. You write a high level configuration file
specifying your inputs and analysis parameters. This input drives a parallel
pipeline that handles distributed execution, idempotent processing restarts and
safe transactional steps. The goal is to provide a shared community resource
that handles the front end data processing component of sequencing analysis,
allowing us to focus on the downstream biology.

## Quick start

1. Install `bcbio-nextgen` and Python dependencies:

          pip install --upgrade bcbio-nextgen
    
2. Edit a [system configuration file][q2] and [sample configuration file][q1]

3. Run analysis, distributed across 8 local cores:

          bcbio_nextgen.py bcbio_system.yaml bcbio_sample.yaml -n 8

[q1]: https://github.com/chapmanb/bcbb/blob/master/nextgen/config/bcbio_sample.yaml
[q2]: https://github.com/chapmanb/bcbb/blob/master/nextgen/config/bcbio_system.yaml

## Documentation

See the [full documentation at ReadTheDocs][d1].

[d1]: https://bcbio-nextgen.readthedocs.org

## Pipelines

### Variant calling

The pipeline implements the [GATK best practice][v1] guidelines for variant
calling, which includes:

- Alignment
- Base Quality Recalibration
- Realignment around indels
- Variant calling. The pipeline supports:
    - [GATK Unified Genotyper][v2] (part of GATK-lite in GATK 2.x)
    - [GATK Haplotype caller][v3] (part of the full, non-open source GATK 2.x)
    - [FreeBayes][v4]
    - [samtools mpileup][v5]
    - [cortex_var][v6]
- Quality filtering, using both [GATK's Variant Quality Score Recalibrator][v7]
  and hard filtering.
- Annotation of effects, using [snpEff][v8]

[v1]: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
[v2]: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
[v3]: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
[v4]: https://github.com/ekg/freebayes
[v5]: http://samtools.sourceforge.net/mpileup.shtml
[v6]: http://cortexassembler.sourceforge.net/index_cortex_var.html
[v7]: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
[v8]: http://snpeff.sourceforge.net/

## Features

### Distributed

The pipeline runs on single multicore machines, in compute clusters managed by
LSF or SGE using [IPython parallel][o8], or on the Amazon cloud.
[This tutorial][o5] describes running the pipeline on Amazon with
[CloudBioLinux][o6] and [CloudMan][o7].

### Galaxy integration

The scripts can be tightly integrated with the [Galaxy][o1]
web-based analysis tool. Tracking of samples occurs via a web based LIMS
system, and processed results are uploading into Galaxy Data Libraries for
researcher access and additional analysis. See the
[installation instructions for the front end][o2] and a
[detailed description of the full system][o3].

[o1]: http://galaxy.psu.edu/
[o2]: https://bitbucket.org/galaxy/galaxy-central/wiki/LIMS/nglims
[o3]: http://bcbio.wordpress.com/2011/01/11/next-generation-sequencing-information-management-and-analysis-system-for-galaxy/
[o5]: http://bcbio.wordpress.com/2011/08/19/distributed-exome-analysis-pipeline-with-cloudbiolinux-and-cloudman/
[o6]: http://cloudbiolinux.org
[o7]: http://wiki.g2.bx.psu.edu/Admin/Cloud
[o8]: http://ipython.org/ipython-doc/dev/index.html
