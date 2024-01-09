REFERENCE='/cluster/projects/pacbio_gsc_gcooper/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
PRIMARY_CONTIGS='/cluster/projects/pacbio_gsc_gcooper/hg38_contigs.bed.gz'
SCRIPT_DIRECTORY='/cluster/home/jlawlor/jasmine_consensus/scripts'
PIPELINE_DIRECTORY='/cluster/home/jlawlor/jasmine_consensus/pipeline'

# Define information about the callers we're putting into jasmine.
# Jasmine sometimes drops INFO descriptions in the headers from the input, so for each caller we 
# make a list with bcftools view -h example.vcf | grep INFO | grep -v '#CHROM' 
ORDERED_CALLERS = [
    {
        'name': 'pbsv',
        'order': 0,
        'info_header': f'{SCRIPT_DIRECTORY}/pbsv_info.txt'
        
    },
    {
        'name': 'sniffles2',
        'order': 1,
        'info_header': f'{SCRIPT_DIRECTORY}/sniffles_info.txt'
        
    }
]
# annotations is a list of INFO fields to transfer or '*' to transfer all INFO fields
# Consider: do I want this to be an input so that all outputs are force regeenrated if we add something new
ANNOTATIONS = [
    {
        'description': 'inhouse_pbsv',
        'vcf': '/cluster/projects/pacbio_gsc_gcooper/resources/pbsv_frequency_20230627/pbsv_single_call_merge_266_individuals_2023-06-27.uniqueid.vcf',
        'annotations': '*'
    },
    # {
    #     'description': 'gnomad4',
    #     'vcf': '/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/gnomad.v4.0.sv.ALL.vcf', # has unique IDs already
    #     'annotations': ['POPMAX_AF','AC','AN','AF']
    # },
    {
        'description': 'hgsvc2',
        'vcf': '/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/hgsvc2_tagged_variants_freeze4_sv_insdel_alt.uniqueid.vcf',
        'annotations': ['AC', 'AN', 'AF', 'SAMPLE']
    },
    {
        'description': 'hprc_giab_pbsv',
        'vcf': '/cluster/projects/pacbio_gsc_gcooper/resources/sv_annotations/original/HPRC_GIAB.GRCh38.pbsv.uniqueid.vcf',
        'annotations': '*'
    }
]

rule all:
    input: 
        expand(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.squashed.vcf', sample=config['samples'].split(',')),
        #expand(PIPELINE_DIRECTORY + '/consensus/{sample}.db/data.mdb', sample=config['samples'].split(',')), 
        #PIPELINE_DIRECTORY + '/annotation.db/data.mdb'

rule prep_pbsv:
    input: '/cluster/projects/pacbio_gsc_gcooper/{sample}/{sample}.sv.vcf.gz'
    output: PIPELINE_DIRECTORY + '/prepped_pbsv/{sample}.pbsv.vcf'
    threads: 1
    resources:
        mem_mb=1*1024
    conda: 'htslib.yaml'
    params:
        primary_contigs=PRIMARY_CONTIGS
    shell:
        '''
        bcftools view -R {params.primary_contigs} -O v -o {output} {input}
        '''

rule prep_sniffles:
    input: '/cluster/projects/pacbio_gsc_gcooper/{sample}/{sample}.snf'
    output: 
        temp_vcf=temp(PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.temp.vcf.gz'),
        vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf'
    threads: 1
    params:
        genome=REFERENCE,
        primary_contigs=PRIMARY_CONTIGS
    resources:
        mem_mb=4*1024
    conda: 'sniffles.yaml'
    shell:
        '''
        sniffles --input {input} --vcf {output.temp_vcf} --reference {params.genome}
        bcftools view -R {params.primary_contigs} -O v -o {output.vcf} {output.temp_vcf}
        '''

rule jasmine_consensus:
    input: 
        pbsv_vcf=PIPELINE_DIRECTORY + '/prepped_pbsv/{sample}.pbsv.vcf',
        sniffles_vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf'
    output: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf',
        #tbi=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz.tbi',
        tmp_vcf=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.tmp.vcf'),
        tmp_dir=temp(directory(PIPELINE_DIRECTORY + '/consensus/{sample}.working/')),
        headers=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.headers.txt')
    threads: 4
    params:
        genome=REFERENCE,
        min_overlap=0.65,
        min_seq_id=0.65,
        headers=' '.join([x['info_header'] for x in ORDERED_CALLERS])
    resources:
        mem_mb=32*1024
    conda: 'jasmine.yaml'
    shell:
        '''
        mkdir -p {output.tmp_dir}
        cat {params.headers} > {output.headers}
        jasmine --min_overlap={params.min_overlap} --min_seq_id={params.min_seq_id} --dup_to_ins --comma_filelist --output_genotypes out_file={output.tmp_vcf} genome_file={params.genome} threads={threads} out_dir={output.tmp_dir} file_list={input.pbsv_vcf},{input.sniffles_vcf}
        bcftools annotate --header-lines {output.headers} {output.tmp_vcf} | bcftools sort -o {output.vcf} 
        '''

rule squash_genotypes:
    input: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf',
    output: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.squashed.vcf',
    conda: 'cyvcf2-lmdbm.yaml' # because it has numpy
    threads: 1
    resources:
        mem_mb=2*1024
    script: 'simple_squash.py'

rule jasmine_annotation_merge:
    input: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz',
        tbi=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz'
    output:
        vcf=PIPELINE_DIRECTORY + '/annotate/{sample}.annomerge.vcf',
        tmp_dir=temp(directory(PIPELINE_DIRECTORY + '/annotate/{sample}.working/'))
    params:
        genome=REFERENCE,
        annotation_vcfs = ','.join([x['vcf'] for x in ANNOTATIONS]),
        min_overlap=0.65,
        min_seq_id=0.65
    threads: 4
    resources:
        mem_mb=48*1024
    conda: 'jasmine.yaml'
    # seems unhappy with duptoins with these inputs
    shell:
        '''
        # these headers will only have the primary source VCF header lines, so they will be complete at this point 
        mkdir -p {output.tmp_dir}
        gunzip -c {input.vcf} > {output.tmp_dir}/input.vcf
        jasmine --require_first_sample --centroid-merging --min_overlap={params.min_overlap} --min_seq_id={params.min_seq_id} --comma_filelist out_file={output.vcf} genome_file={params.genome} threads={threads} out_dir={output.tmp_dir} file_list={output.tmp_dir}/input.vcf,{params.annotation_vcfs}
        '''

rule annotate_jasmine:
    input:
        vcf=PIPELINE_DIRECTORY + '/annotate/{sample}.annomerge.vcf',
        annotation_db=PIPELINE_DIRECTORY + '/annotation.db/data.mdb'
    output:
        vcf=PIPELINE_DIRECTORY + '/annotate/{sample}.annomerge.annotated.vcf'
    params:
        annotations = ANNOTATIONS
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=8*1024
    script: 'transfer_annotations.py'

# This takes about 15 minutes for inhouse + gnomad + hgsvc2 + hprc_giab
# TODO: Figure out how to make this always rebuild if the parameters change
rule build_annotation_table:
    output: PIPELINE_DIRECTORY + '/annotation.db/data.mdb'
    params:
        annotations = ANNOTATIONS
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=4*1024
    script: 'build_annotation_table.py'


rule sample_annotation_table:
    input:
        sniffles_vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf',
        pbsv_vcf=PIPELINE_DIRECTORY + '/prepped_pbsv/{sample}.pbsv.vcf'
    output: PIPELINE_DIRECTORY + '/consensus/{sample}.db/data.mdb'
    params:
        annotations = [
            {
                'description': 'sniffles',
                'vcf': 'sniffles_vcf',
                'annotations': '*'
            },
            {
                'description': 'pbsv',
                'vcf': 'pbsv_vcf',
                'annotations': '*'
            }
        ]
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=16*1024
    script: 'build_annotation_table.py'