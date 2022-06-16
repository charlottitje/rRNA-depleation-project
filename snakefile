SETNAME= config["set"]

SPECIES= glob_wildcards("genomes/{organism}.fna").organism


rule ALL:
    input: 
        expand("{set_name}.order.txt", set_name= SETNAME),
        expand("{set_name}.tsv", set_name= SETNAME)
      
rule barrnapFASTA:
    input:
        "genomes/{species}.fna"
    output:
        fasta="rRNA_fasta/{species}.fasta"
    log:
        err="logs/barnapFASTA.{species}.err"
    shell:
        "barrnap --kingdom euk --quiet --outseq {output.fasta} < {input} 2> {log.err}"
        
rule barrnapGFF:
    input:
        "genomes/{species}.fna"
    output:
        gff="rRNA_gff/{species}.gff"
    log:
        err="logs/barnapGFF.{species}.err"
    shell:
        "barrnap --kingdom euk --quiet  < {input} >{output.gff} 2> {log.err}"
        
rule guide_finder:
    input:
        "rRNA_fasta/{species}.fasta"
    output:
        "guide_finder/{species}.fasta"
    log:
        out="logs/guide_finder.{species}.out",
        err="logs/guide_finder.{species}.err"
    shell:
        "python3 scripts/guide_finder.py -i {input} -o {output}"
        
rule bowtie_index:
    input:
        genome="genomes/{index}.fna"
    output:
        touch("BT_index/touch_files/{index}")
    params:
        idx="BT_index/{index}"        
    log:
        out="logs/bowtie-index.{index}.out",
        err="logs/bowtie-index.{index}.err"
    shell:
        "bowtie-build {input.genome} {params.idx}"

rule bowtie_aling:
    input:
        guides="guide_finder/{species}.fasta",
        idx="BT_index/touch_files/{index}"
    output:
        "bowtie_out/{index}.{species}.sam"
    params:
        ref_idx= "BT_index/{index}"
    log:
        out="logs/bowtie-aling.{index}.{species}.out",
        err="logs/bowtie-aling.{index}.{species}.err"
    shell:
        "bowtie -a -S -v 3 -f {params.ref_idx} {input.guides} {output} > {log.out} 2> {log.err}"
        
rule guide_curator:
    input:
        gff= expand("rRNA_gff/{species}.gff", species=SPECIES),
        sam= expand("bowtie_out/{index}.{species}.sam", species= SPECIES, index= SPECIES),
        fasta= expand("guide_finder/{species}.fasta", species=SPECIES)
    output:
        order= "{set_name}.order.txt",
        tsv= "{set_name}.tsv"
    params:
        gffdir= "rRNA_gff/",
        samdir= "bowtie_out/"
    log:
        out="logs/bowtie-index.{set_name}.out",
        err="logs/bowtie-index.{set_name}.err"
    shell:
        "python3 scripts/guide_curator.py --gff {params.gffdir} --sam {params.samdir} --set {wildcards.set_name}"
