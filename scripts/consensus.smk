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
        'order': 1,
        'info_header': f'{SCRIPT_DIRECTORY}/pbsv_info.txt'
        
    },
    {
        'name': 'sniffles2',
        'order': 2,
        'info_header': f'{SCRIPT_DIRECTORY}/sniffles_info.txt'
        
    }
]

rule consensus:
    input: expand(PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz', sample=config['samples'].split(','))

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
        vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz',
        tbi=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf.gz.tbi',
        tmp_vcf=PIPELINE_DIRECTORY + '/consensus/{sample}.consensus.sv.vcf',
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
        jasmine --leave_breakpoints --min_overlap={params.min_overlap} --min_seq_id={params.min_seq_id} --dup_to_ins --comma_filelist --output_genotypes out_file={output.tmp_vcf} genome_file={params.genome} threads={threads} out_dir={output.tmp_dir} file_list={input.pbsv_vcf},{input.sniffles_vcf}
        bcftools annotate --header-lines {output.headers} {output.tmp_vcf} | bcftools sort -o {output.vcf} 
        tabix -p vcf {output.vcf}
        '''