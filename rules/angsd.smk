rule angsd_GL2:
  input:
    ref = config["ref_rapid"],
    bamlist = 'depth/stats/realignedBAM_df3.list'
  output:
    touch('angsd/angsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/angsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Calculate genotype likelihoods and call SNPs for all bam files using angsd setting various cutoffs """
  shell:
    """
    module load angsd/0.935
    angsd -b {input.bamlist} -ref {input.ref} -out angsd/angsd_GL2_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 30 -minQ 20 -skipTriallelic 1 -minInd {wildcards.IND} -setMinDepth {wildcards.MinDepth} -setMaxDepth {wildcards.MaxDepth} -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf {wildcards.minMaf} -doCounts 1 -nThreads {threads} 2> {log}
    """
