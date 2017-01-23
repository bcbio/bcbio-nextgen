NAMES = {
    'lane': 'Test1',
    'lb': None,
    'pu': 'Test1',
    'sample': 'Test1',
    'rg': 'Test1',
    'pl': 'illumina'
}

DATA = {
    'files': [
        '/bcbio-nextgen/tests/test_automated_output/trimmed/1_1_Test1.trimmed.fq.gz',
        '/bcbio-nextgen/tests/test_automated_output/trimmed/1_2_Test1.trimmed.fq.gz'
    ],
    'dirs': {
        'config': '/bcbio-nextgen/tests/test_automated_output',
        'fastq': '/bcbio-nextgen/tests/data/test_fusion',
        'work': '/bcbio-nextgen/tests/test_automated_output',
        'flowcell': '/bcbio-nextgen/tests/data/test_fusion',
        'galaxy': '/bcbio-nextgen/tests/data/automated'
    },
    'lane': '1',
    'description': 'Test1',
    'reference': {
        'genome_context': [
            '/bcbio-nextgen/tests/data/genomes/hg19/coverage/problem_regions/GA4GH/test.bed.gz',
            '/bcbio-nextgen/tests/data/genomes/hg19/coverage/problem_regions/GA4GH/test2.bed.gz'
        ],
        'fasta': {
            'base': '/bcbio-nextgen/tests/data/genomes/hg19/seq/hg19.fa'
        },
        'star': {
            'indexes': [
                '/bcbio-nextgen/tests/data/genomes/hg19/star/chrLength.txt',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/sjdbList.out.tab',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/SA',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/Genome',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/SAindex',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/chrStart.txt',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/chrName.txt',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/chrNameLength.txt',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/genomeParameters.txt',
                '/bcbio-nextgen/tests/data/genomes/hg19/star/Log.out'
            ]
        },
        'rtg': '/bcbio-nextgen/tests/data/genomes/hg19/rtg/hg19.sdf'
    },
    'sam_ref': '/bcbio-nextgen/tests/data/genomes/hg19/seq/hg19.fa',
    'genome_resources': {
        'rnaseq': {
            'transcripts': '/bcbio-nextgen/tests/data/genomes/hg19/rnaseq/ref-transcripts.gtf',
            'transcripts_mask': '/bcbio-nextgen/tests/data/genomes/hg19/rnaseq/ref-transcripts-mask.gtf',
            'gene_bed': '/bcbio-nextgen/tests/data/genomes/hg19/rnaseq/ref-transcripts.bed'
        },
        'version': 7,
        'variation': {
            'train_omni': '/bcbio-nextgen/tests/data/genomes/hg19/variation/1000G_omni2.5.vcf.gz',
            'dbnsfp': '/bcbio-nextgen/tests/data/genomes/hg19/variation/dbNSFP_v2.5.gz',
            'cosmic': '/bcbio-nextgen/tests/data/genomes/hg19/variation/cosmic-v68-hg19.vcf.gz',
            'ancestral': '/bcbio-nextgen/tests/data/genomes/hg19/variation/human_ancestor.fa',
            'train_hapmap': '/bcbio-nextgen/tests/data/genomes/hg19/variation/hapmap_3.3.vcf.gz',
            'train_1000g': '/bcbio-nextgen/tests/data/genomes/hg19/variation/1000G_phase1.snps.high_confidence.vcf.gz',
            'dbsnp': '/bcbio-nextgen/tests/data/genomes/hg19/variation/dbsnp_132.vcf.gz',
            'train_indels': '/bcbio-nextgen/tests/data/genomes/hg19/variation/Mills_Devine_2hit.indels.vcf.gz'
        },
        'srnaseq': {
            'srna-trasncripts': '/bcbio-nextgen/tests/data/genomes/hg19/srnaseq/srna-transcripts.gtf',
            'mirbase': '/bcbio-nextgen/tests/data/genomes/hg19/srnaseq/hairpin.fa'
        },
        'aliases': {
            'snpeff': 'hg19',
            'human': True,
            'ensembl': 'homo_sapiens_vep_83_GRCh37'
        }
    },
    'provenance': {
        'data': '/bcbio-nextgen/tests/test_automated_output/provenance/data_versions.csv',
        'entity': 'bcdd2c84-b800-11e6-a323-0242ac110002.prepare_sample.0.trim_sample.0.process_alignment.0',
        'db': None,
        'programs': '/bcbio-nextgen/tests/test_automated_output/provenance/programs.txt'
    },
    'rgnames': {
        'lane': 'Test1',
        'lb': None,
        'pu': 'Test1',
        'sample': 'Test1',
        'rg': 'Test1',
        'pl': 'illumina'
    },
    'upload': {
        'dir': '/bcbio-nextgen/tests/test_automated_output/upload',
        'run_id': ''
    },
    'analysis': 'RNA-seq',
    'name': ['', 'Test1'],
    'genome_build': 'hg19',
    'config': {
        'galaxy_config': '/bcbio-nextgen/tests/data/automated/universe_wsgi.ini',
        'resources': {
            'gatk': {
                'jvm_opts': ['-Xms500m', '-Xmx3500m']
            },
            'default': {
                'cores': 16,
                'jvm_opts': ['-Xms750m', '-Xmx3500m'],
                'memory': '3G'
            },
            'express': {'memory': '8g'},
            'seqcluster': {'memory': '8g'},
            'program_versions': '/bcbio-nextgen/tests/test_automated_output/provenance/programs.txt',
            'dexseq': {'memory': '10g'},
            'macs2': {'memory': '8g'},
            'snpeff': {
                'jvm_opts': ['-Xms750m', '-Xmx4g']
            },
            'qualimap': {'memory': '4g'}
        },
        'log_dir': '/var/log/bcbio',
        'algorithm': {
            'nomap_split_targets': 200,
            'trim_reads': 'read_through',
            'qc': ['fastqc', 'qualimap_rnaseq', 'samtools', 'gemini'],
            'archive': [],
            'recalibrate': False,
            'mark_duplicates': True,
            'nomap_split_size': 250,
            'quality_format': 'illumina',
            'aligner': 'star',
            'validate_regions': None,
            'realign': False,
            'tools_off': [],
            'fusion_mode': True,
            'variant_regions': None,
            'coverage_interval': None,
            'adapters': ['truseq', 'polya'],
            'validate': None, 'num_cores': 1, 'tools_on': []
        },
        'bcbio_system': '/bcbio-nextgen/tests/test_automated_output/bcbio_system-merged.yaml'
    },
    'resources': {},
    'metadata': {
        'batch': None,
        'phenotype': ''
    }
}

