rule raxml:
    """

    Creates a tree using RAxML

    """
    input:
        aln = "output/align/macse/concat/{clade}.aligned.cat.fasta"
    output:
        outdir=directory("output/phylo/raxml/raxml_{clade}"),
        bootstraps="output/phylo/raxml/raxml_{clade}/RAxML_bootstrap.{clade}_out",
        tree="output/phylo/raxml/raxml_{clade}/RAxML_bestTree.{clade}_out"
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
        "output/logs/phylo/raxml/{clade}-raxml.log"
    resources:
        mem_mb=config['mem_mb']['raxml']
    shell:
        """
        function realpath {{ echo $(cd $(dirname $1); pwd)/$(basename $1); }}
        outdir=`realpath {output.outdir}`
        mkdir -p $outdir
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
        {params.other} 2> {log} 1>&2
        """
