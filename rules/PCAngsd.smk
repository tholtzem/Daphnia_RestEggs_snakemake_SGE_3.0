
rule PCAngsd_covmat:
  input: 
    touched = 'angsd/angsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    touch('pcangsd/PCAngsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/PCAngsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_GL2_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz -o pcangsd/PCAngsd_GL2_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_covmat 2> {log}
    """

rule plot_covMat:
  input:
    #touched = 'pcangsd/PCAngsd_GL2_{IND}.done',
    #covMat ='pcangsd/PCAngsd_GL2_{IND}_covmat.cov'
    touched = 'pcangsd/PCAngsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done',
    metadata = 'list/samples184_metadata.tsv'
  output:
    pdf = 'pcangsd/PCAngsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_covmat.pdf'
  log: 'log/PCAngsd_plotcovmat_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat.R pcangsd/PCAngsd_GL2_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_covmat.cov {input.metadata} {output.pdf} 2> {log}
    """

