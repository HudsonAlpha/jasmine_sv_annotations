import json 

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
with open(SCRIPT_DIRECTORY + '/annotations.json') as f:
    ANNOTATIONS = json.load(f)

rule all:
    input: 
        expand(PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.sv.vcf.gz', sample=config['samples'].split(',')),
        #PIPELINE_DIRECTORY + 'ANNOTATIONS.json',
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
        sniffles_vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf',
    output: 
        vcf=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.unsquashed.vcf'),
        tmp_vcf=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.tmp.vcf'),
        tmp_dir=temp(directory(PIPELINE_DIRECTORY + '/consensus/{sample}.working/')),
        headers=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.headers.txt'),
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

rule annotate_consensus:
    input: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.unsquashed.vcf',
        pbsv_vcf=PIPELINE_DIRECTORY + '/prepped_pbsv/{sample}.pbsv.vcf',
        sniffles_vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf',
        sample_annotation_table=PIPELINE_DIRECTORY + '/consensus/{sample}.db/data.mdb',
    output:
        vcf=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.unsquashed.annotated.vcf')
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=8*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/transfer_annotations.py --input {input.vcf} --output {output.vcf} --annotation_db {input.sample_annotation_table} --vcftuples {input.pbsv_vcf},pbsv {input.sniffles_vcf},sniffles
        '''

rule squash_genotypes:
    input: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.unsquashed.annotated.vcf'
    output: 
        vcf=temp(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf'),
    conda: 'cyvcf2-lmdbm.yaml' # because it has numpy
    threads: 1
    resources:
        mem_mb=2*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/simple_squash.py --input {input.vcf} --output {output.vcf} --sample {wildcards.sample}
        '''

rule clean_bnds:
    input:
        vcf=PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.vcf',
    output:
        vcf=temp(PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.cleanbnd.vcf')
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=2*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/clean_bnds.py --input {input.vcf} --output {output.vcf}
        '''

rule jasmine_annotation_merge:
    input: 
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf',
    output:
        vcf=PIPELINE_DIRECTORY + '/annotated/{sample}.annomerge.vcf',
        tmp_dir=temp(directory(PIPELINE_DIRECTORY + '/annotated/{sample}.working/'))
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
        jasmine --require_first_sample --centroid-merging --min_overlap={params.min_overlap} --min_seq_id={params.min_seq_id} --comma_filelist out_file={output.vcf} genome_file={params.genome} threads={threads} out_dir={output.tmp_dir} file_list={input.vcf},{params.annotation_vcfs}
        '''

rule transfer_consensus_annotations:
    input:
        vcf=PIPELINE_DIRECTORY + '/annotated/{sample}.annomerge.vcf',
        annotation_db=PIPELINE_DIRECTORY + '/annotation.db/data.mdb',
    params:
        annotation_json=SCRIPT_DIRECTORY + '/annotations.json'
    output:
        vcf=temp(PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.vcf')
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=8*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/transfer_annotations.py --input {input.vcf} --output {output.vcf} --annotation_db {input.annotation_db} --annotation_json {params.annotation_json}
        '''
# This takes about 15 minutes for inhouse + gnomad + hgsvc2 + hprc_giab

rule bgzip_final:
    input: 
        vcf=PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.cleanbnd.vcf'
    output: 
        vcf=PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.sv.vcf.gz',
        tbi=PIPELINE_DIRECTORY + '/annotated/{sample}.consensus.sv.vcf.gz.tbi'
    conda: 'htslib.yaml'
    threads: 1
    resources:
        mem_mb=1*1024
    shell:
        '''
        bcftools sort -O z -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        '''

rule build_annotation_table:
    input: SCRIPT_DIRECTORY + '/annotations.json'
    output: PIPELINE_DIRECTORY + '/annotation.db/data.mdb'
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=24*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/build_annotation_table.py --output {output} --annotations {input}
        '''


rule sample_annotation_table:
    input:
        sniffles_vcf=PIPELINE_DIRECTORY + '/prepped_sniffles/{sample}.sniffles.vcf',
        pbsv_vcf=PIPELINE_DIRECTORY + '/prepped_pbsv/{sample}.pbsv.vcf',
    output: 
        anno_table=PIPELINE_DIRECTORY + '/consensus/{sample}.db/data.mdb',
    conda: 'cyvcf2-lmdbm.yaml'
    threads: 1
    resources:
        mem_mb=16*1024
    shell:
        '''
        python {SCRIPT_DIRECTORY}/build_annotation_table.py --vcftuples {input.pbsv_vcf},pbsv {input.sniffles_vcf},sniffles --output {output.anno_table}
        '''