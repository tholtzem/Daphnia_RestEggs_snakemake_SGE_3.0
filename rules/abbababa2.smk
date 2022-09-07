rule prepare_bamlist:
  input:
    pop1 = 'list/pop_list/{POP1}.txt',
    pop2 = 'list/pop_list/{POP2}.txt',
    pop3 = 'list/pop_list/{POP3}.txt',
    pop4 = 'list/pop_list/{POP4}.txt'
  output:
    #bamlist = 'list/abbababa/bam.list',
    #popsize = 'list/abbababa/sizefile.list'
    pop1 = touch('list/abbababa/{POP1}.added2bamlist.done'),
    pop2 = touch('list/abbababa/{POP2}.added2bamlist.done'),
    pop3 = touch('list/abbababa/{POP3}.added2bamlist.done'),
    pop4 = touch('list/abbababa/{POP4}.added2bamlist.done')
  log:
    'log/abbababa/added2bamlist.log'
  threads: 12
  message:
    """ Preparea list of bam files required for the ABBA-BABA test (D-statistic) """
  shell:
    """
    POP1=$(echo {input.pop1} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP2=$(echo {input.pop2} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP3=$(echo {input.pop3} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP4=$(echo {input.pop4} | sed 's!df_!!g' | sed 's!.txt!!g)

    cat {input.pop1} > list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop2} >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop3} >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop4} >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list 2> {log}
    """

rule prepare_sizefile:
  input:
    pop1 = 'list/pop_list/{POP1}.txt',
    pop2 = 'list/pop_list/{POP2}.txt',
    pop3 = 'list/pop_list/{POP3}.txt',
    pop4 = 'list/pop_list/{POP4}.txt'
  output:
    #bamlist = 'list/abbababa/bam.list',
    #popsize = 'list/abbababa/sizefile.list'
    pop1 = touch('list/abbababa/{POP1}.added2sizefile.done'),
    pop2 = touch('list/abbababa/{POP2}.added2sizefile.done'),
    pop3 = touch('list/abbababa/{POP3}.added2sizefile.done'),
    pop4 = touch('list/abbababa/{POP4}.added2sizefile.done')
  log:
    'log/abbababa/added2sizefile.log'
  threads: 12
  message:
    """ Preparea list of bam files required for the ABBA-BABA test (D-statistic) """
  shell:
    """
    POP1=$(echo {input.pop1} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP2=$(echo {input.pop2} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP3=$(echo {input.pop3} | sed 's!df_!!g' | sed 's!.txt!!g)
    POP4=$(echo {input.pop4} | sed 's!df_!!g' | sed 's!.txt!!g)

    cat {input.pop1} | wc -l > list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop2} | wc -l >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop3} | wc -l >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list &&
    cat {input.pop4} | wc -l >> list/abbababa/${POP1}_${POP2}_${POP3}_${POP4}_bam.list 2> {log}
    """


rule ABBA_BABA:
  input:
    ref = config["ref_rapid"],
    bamlist = 'list/abbababa/bam.list'
    popsize = 'list/abbababa/sizefile.list'
   'list/abbababa/{POP1}.added2bamlist.done'
   'list/abbababa/{POP2}.added2bamlist.done'
   'list/abbababa/{POP3}.added2bamlist.done'
   'list/abbababa/{POP4}.added2bamlist.done'
   'list/abbababa/{POP1}.added2sizefile.done'
   'list/abbababa/{POP2}.added2sizefile.done'
   'list/abbababa/{POP3}.added2sizefile.done'
   'list/abbababa/{POP4}.added2sizefile.done'
  output:
   touch('abbababa/abbababa.done')
  log:
    'log/abbababa/abbababa.log'
  threads: 12
  message:
    """ Compute ABBA-BABA test (D-statistic) using angsd, allows for multiple individuals in each group """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -doAbbababa2 -b {input.bamlist} -useLast 1 -sizeFile {input.popsize} -out abbababa/abbababa -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -doCounts 1 2> {log}
    """
