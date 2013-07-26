Pipelines
---------

Variant calling
~~~~~~~~~~~~~~~

bcbio-nextgen implements configurable best-practice pipelines for SNP
and small indel calling:

-  Sequence alignment:

   - `bwa mem`_
   - `novoalign`_
   - `bowtie2`_
   - `mosaik`_

-  Base Quality Recalibration
-  Realignment around indels
-  Variant calling:

   -  `GATK Unified Genotyper`_ (supports both GATK-lite in GATK 2.3
      and commercial restricted version in GATK 2.4+)
   -  `GATK Haplotype caller`_ (part of the commercially restricted GATK 2.4+)
   -  `FreeBayes`_
   -  `samtools mpileup`_
   -  `cortex\_var`_

-  Quality filtering, using either
   `GATK's Variant Quality Score Recalibrator`_ or hard filtering.
-  Annotation of variant effects, using `snpEff`_
-  Variant exploration and prioritization, using `GEMINI`_

It follows approaches from:

- `GATK best practice`_ guidelines for variant calling
- Marth Lab's `gkno pipelines`_

.. _GATK best practice: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
.. _GATK Unified Genotyper: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
.. _GATK Haplotype caller: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
.. _FreeBayes: https://github.com/ekg/freebayes
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _cortex\_var: http://cortexassembler.sourceforge.net/index_cortex_var.html
.. _GATK's Variant Quality Score Recalibrator: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
.. _snpEff: http://snpeff.sourceforge.net/
.. _bwa mem: http://bio-bwa.sourceforge.net/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _mosaik: https://github.com/wanpinglee/MOSAIK
.. _novoalign: http://www.novocraft.com
.. _gkno pipelines: http://gkno.me/pipelines.html
.. _GEMINI: http://gemini.readthedocs.org/en/latest/

