rule bb_refIndex:
  input:
    ref = config['ref_HiC']
  output:
    touch('ref/bb_indexRef.done')
  log: 'log/bb_indexRef.log'
  threads: 2
  message: """ --- Index reference genome for bbmap --- """
  shell:
    """
    bbmap.sh -Xmx4g t={threads} ref={input.ref} 2> {log}
    """



rule bbmap:
  input:
    kraken_r1 = lambda wildcards: getKrakenHome(wildcards.sample)[0],
    kraken_r2 = lambda wildcards: getKrakenHome(wildcards.sample)[1],
    ref = config['ref_HiC'],
    idx = 'ref/bb_indexRef.done',
    SAMPLETABLE = 'list/eggsLC_pelagial2021_metadata.tsv'
    #'ref/genome/1/summary.txt'
  output:
    bam = 'bbmap/HiC/minid95/{sample}.bam',
    bhist = 'bbmap/HiC/minid95/stats/bhist/{sample}.bhist.txt',
    qhist = 'bbmap/HiC/minid95/qhist/{sample}.qhist.txt',
    lhist = 'bbmap/HiC/minid95/lhist/{sample}.lhist.txt',
    covstats = 'bbmap/HiC/minid95/cov/{sample}.covstats.txt',
    covhist = 'bbmap/HiC/minid95/cov/{sample}.covhist.txt',
    basecov = 'bbmap/HiC/minid95/cov/{sample}.basecov.txt',
    bincov = 'bbmap/HiC/minid95/cov/{sample}.bincov.txt'
  log: 'log/HiC/minid95/bbmap_HiC_{sample}_bam.log'
  threads: 24
  message: """ --- Mapping reads to reference genome with minid 0.95, convert 2 bam, exclude unmapped reads, only keep reads with minq => 20 --- """
  shell:
    """
    id=`echo {input.kraken_r1} | sed -e "s/_R1.trmdfilt.keep.fq.gz$//" | cut -f 2 -d '/'`
    echo $id
    
    RG_ID=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 2`
    RG_LB=LIB1of_${{RG_ID}}
    RG_SM=$RG_ID
    RG_PL=ILLUMINA
    RG_PU=LIB1of_${{RG_ID}}
    
    echo $RG_ID
    echo $RG_LB
    echo $RG_SM
    echo $RG_PL
    echo $RG_PU
    
    /usr/site/hpc/bin/sysconfcpus -n 24 bbmap.sh -Xmx200g t={threads} ref={input.ref} in1={input.kraken_r1} in2={input.kraken_r2} out=stdout.sam minid=0.95 k=13 bw=0 ordered=t rgid=$RG_ID rglb=$RG_LB rgsm=$RG_SM rgpl=$RG_PL rgpu=$RG_PU overwrite=f unpigz=t bhist={output.bhist} qhist={output.qhist} lhist={output.lhist} covstats={output.covstats} covhist={output.covhist} basecov={output.basecov} bincov={output.bincov} | samtools view -F 4 -Shu -q 20 | samtools sort - -o {output.bam} 2> {log}
    """


rule bamIndex:
  input:
    'bbmap/HiC/minid95/{sample}.bam'
  output:
    'bbmap/HiC/minid95/{sample}.bam.bai'
  log: 'log/HiC/minid95/index_HiC_bam{sample}.log'
  threads: 2
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input} {output} 2> {log}
    """


rule remove_duplicates:
  input:
    bam = "bbmap/rapid/{sample}.bam",
    bai = "bbmap/rapid/{sample}.bam.bai"
  output:
    deDup = 'deDup/{sample}.dedup.bam',
    metrics = 'deDup/{sample}.dedup.metrics.txt'
  log: 'log/{sample}.dedup.bam.log'
  threads: 12
  message: """--- Removing duplicates of bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs2.0/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.bam} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
    """


rule samtools_minq20:
  input:
    '/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/deDup/{sample}.dedup.bam'
  output:
    'deDup/{sample}.dedup.bam'
  log: 'log/{sample}.ref_minq20.log'
  message: """ --- Get deduplicated reference clones, exclude unmapped reads and only keep reads with minq => 20 --- """
  shell:
    """
    samtools view -F 4 -Shu -q 20 {input} | samtools sort - -o {output} 2> {log}
    """


rule bamIndex_Refclones:
  input:
    bam = 'deDup/{sample}.dedup.bam',
    logs = 'log/{sample}.ref_minq20.log'
  output:
    'deDup/{sample}.dedup.bam.bai'
  threads: 2
  log: 'log/{sample}.dedup.Indexref.log'
  message: """ --- Indexing reference clones with samtools --- """
  shell:
    """
    samtools index {input.bam} {output} 2> {log}
    """


rule clip_overlap:
  input:
    deDup = 'deDup/{sample}.dedup.bam'
    #log_files = 'log/{sample}.dedup.Indexref.log'
  output:
    clip = 'deDup/{sample}.overlapclipped.bam' 
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 12
  message:
    """ Clip overlapping paired end reads """
  shell:
    """
    bam clipOverlap --in {input.deDup} --out {output.clip} --stats 2> {log}
    """


rule Index_clippedBAM:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    idx = 'deDup/{sample}.overlapclipped.bam.bai'
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 2
  message: """--- Indexing clipped BAM files with samtools ---"""
  shell:
    """
    samtools index {input.clip} {output.idx} 2> {log}
    """


rule refIndex:
  input:
    ref = config['ref_rapid']
  output:
    "ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai"
  log: 'log/dgal_ra_refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """


rule ref_Dict:
  input:
    ref = config['ref_rapid']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  log: 'log/dgal_ra_refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs2.0/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """

    
rule ls_ClipBam:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    touch('deDup/{sample}.added2ClippedList.done')
  log: 'log/{sample}.added2ClippedList.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls {input.clip} >> list/overlapclippedBAM.list 2> {log}
    """


rule list_indels:
  input:
    clip = 'list/overlapclippedBAM.list',
    ref = config['ref_rapid'],
    idx_ref = "ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai",
    dict_ref = 'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  output:
    indels = 'list/indels.list'
  log: 'log/listIndels.log'
  threads: 48
  message:
    """ Create list of potential indels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs2.0/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 48 java -Xmx440g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output.indels} -drf BadMate 2> {log}
    """


rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config['ref_rapid'],
    indels = 'list/indels.list'
  output:
    realigned = 'realigned/{sample}.realigned.bam'
  log: 'log/{sample}.realigned.bam.log'
  threads: 12
  message:
    """ Realign in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs2.0/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.indels} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """


rule samtools_depth:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    depth = 'depth/{sample}.realigned.bam.depth.gz'
  log: 'log/{sample}.realigned.bam.depth.log'
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth} 2> {log}
    """


rule ls_depth:
  input:
    'depth/{sample}.realigned.bam.depth.gz'
  output:
    touch('depth/{sample}.added2DepthList.list')
  log: 'log/{sample}.depth.list.log'
  message: """ --- Creating a sample list of samtools-depth files for Rscript --- """
  shell:
    """
    ls {input} | cut -f2 -d '/' >> depth/depth.list 2> {log}
    """


rule read_depth:
  input:
    'depth/depth.list'
  output:
    'depth/stats/depth_statistics.txt'
    #touch('depth/stats/genome_stats.done')
  log: 'log/genome_stats.log'
  threads: 12
  message:
    """ --- Running Rscript to plot the genome-wide distribution of coverage --- """
  shell:
    """
    Rscript scripts/read_depth.R {input} {output} 2> {log} 
    """


rule plot_summary:
  input:
    'depth/stats/depth_statistics_{sets}.txt'
  output:
    touch('depth/plots/plot_summary_{sets}.done')
  log: 'log/plot_summary_{sets}.log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/plot_summary.R {input} {wildcards.sets} 2> {log}
    """


rule rbind_depthFilter:
  input:
    args1 = 'depth/stats/LC_depthFilter.list',
    args2 = 'depth/stats/ALL_depthFilter.list',
    args3 = 'depth/stats/LC_REF_depthFilter.list',
    args4 = 'depth/stats/LC_withoutREF_depthFilter.list'
  output:
    'depth/stats/depthFilter.list'
  log: 'log/rbind_dfFilter.log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/rbind_dfFilter.R {input.args1} {input.args2} {input.args3} {input.args4} {output} 2> {log}
    """


rule new_metadata_table:
  input:
    'depth/plots/plot_summary_{sets}.done'#,
    #'depth/stats/realignedBAM_df3.list'
  output:
    'list/new_metadata_table_{sets}.done'
  message: """ --- Create a metadata table for the samples passing the mean depth treshold ( > 1) --- """
  shell:
    """
    ./scripts/create_newmetadata_list.sh 2> {log} 
    """


rule genome_coverage_bed:
  input:
    'realigned/{sample}.realigned.bam'
  output:
    'bedtools/{sample}.realigned.genomecov.bed'
  log: 'log/{sample}.realigned.genomecov.log'
  threads: 12
  message:
    """ Computes BED summaries using bedtools """
  shell:
    """
    bedtools genomecov -ibam {input} > {output} 2> {log}
    """


rule plot_gencov:
  input:
    bed= 'bedtools/{sample}.realigned.genomecov.bed'
  output:
    pdf = 'bedtools/plots/{sample}.realigned.genomecov.pdf'
  log: 'log/{sample}.realigned.genomecov_plot.log'
  threads: 4
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/plot_gene_covs.R {input.bed} {output.pdf} 2> {log}
    """
