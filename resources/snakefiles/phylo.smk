def get_alignment(clade):
    align_fp = clades_df.loc[clade, 'alignment']
    return(align_fp)


rule raxml:
    """

    Creates a tree using RAxML

    """
    input:
        aln = lambda wildcards: get_alignment(wildcards.clade)
    output:
        directory("output/phylo/raxml/raxml_{clade}")
    params:
        raxml=config['raxml']['executable'],
        temp_dir=directory("output/{clade}_temp/"),
        model=config['raxml']['model'],
        algorithm=config['raxml']['algorithm'],
        seed=config['raxml']['seed'],
        bootstraps=config['raxml']['bootstraps'],
        other=config['raxml']['other']
    conda:
        "../env/raxml.yaml"
    threads:
        config['threads']['raxml']
    benchmark:
        "output/benchmarks/phylo/raxml/{clade}-raxml.txt"
    log:
        "output/benchmarks/phylo/raxml/{clade}-raxml.log"
    resources:
        mem_mb=config['mem_mb']['raxml']
    shell:
        """
        outdir=`readlink -f {output}`
        {params.raxml} \
        -s {input.aln} \
        -w $outdir \
        -n {wildcards.clade}_out \
        -m {params.model} \
        -f {params.algorithm} \
        -x {params.seed} \
        -N {params.bootstraps} \
        -p {params.seed} \
        -T {threads} \
        {params.other}
        """
