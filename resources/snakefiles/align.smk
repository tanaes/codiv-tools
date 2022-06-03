import pandas as pd
from glob import glob
from os.path import join, basename

"""
config file: has target genes
sample file: has sample name / path to sample file / clade
gene metadata file: has gene name info

- for a clade:
    - per gene:
        - (in python) grabs all the available genes and makes a fasta
        - (in snakemake) aligns genes
    - concatenates alignments

"""

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def clean_seq(seq):
    full_pattern = re.compile('[^a-zA-Z\-]')

    return re.sub(full_pattern, '-', seq)

def read_fasta_to_series(fp):
    fasta_dict = {}
    with open(fp, 'r') as f:
        for name, seq, _ in readfq(f):
            fasta_dict[name] = seq
    return(fasta_dict)

def load_alignments(aln_files):

    len_dict = {}
    aln_dict = {}
    for aln in aln_files:
        gene_name = basename(aln)[:-9]
        aln_dict[gene_name] = read_fasta_to_series(aln)
        len_dict[gene_name] = len(aln_dict[gene_name][list(aln_dict[gene_name].keys())[0]])
    aln_df = pd.DataFrame.from_dict(aln_dict, orient='index')
    len_series = pd.Series(len_dict)

    return(aln_df, len_series)

def cat_alignments(aln_df, len_series):
    fasta_str = ''
    partitions = ''
    cat_lens = []
    for genome in aln_df.columns:
        cat_seq = ''
        for gene in aln_df.index:
            seq = aln_df.loc[gene, genome]
            if pd.isna(seq):
                seq = '-'*len_series[gene]
            cat_seq += clean_seq(seq)
        fasta_str += '>{0}\n{1}\n'.format(genome, cat_seq)
        cat_lens.append(len(cat_seq))
        
    last = 1
    for gene in len_series.index:
        
        gene_len = len_series[gene]
        if gene_len % 3 != 0:
            print('warning, {0} is not divisible by 3'.format(gene))
        partitions += ('DNA, {0}codon12 = {1}-{4}\\3, {2}-{4}\\3\n'
                       'DNA, {0}codon3 = {3}-{4}\\3\n'.format(gene.replace('.','_'),
                                                             last,
                                                             last+1,
                                                             last+2,
                                                             int(last+gene_len)))
        last += gene_len
    
    aln_lens = pd.Series(cat_lens, index=aln_df.columns)
        
    if len(set(aln_lens)) != 1:
        print('warning, not all alignments are same length')
    
    return(fasta_str, partitions, aln_lens)


rule align_seqs_macse:
    """

    Aligns a set of sequences using MACSE

    """
    input:
        fasta_dir=lambda wildcards: clades_df.loc[str(wildcards.clade), 'fasta_dir']
    output:
        outdir=directory("output/align/macse/{clade}"),
        done=touch("output/align/macse/{clade}/{clade}.done")
    params:
        gc=lambda wildcards: clades_df.loc[str(wildcards.clade),'genetic_code'],
        fasta_ext=config['macse']['fasta_ext'],
        other=config['macse']['other'],
    conda:
        "../env/align.yaml"
    threads:
        config['threads']['macse']
    benchmark:
        "output/benchmarks/align/macse/{clade}-macse.txt"
    log:
        "output/logs/align/macse/{clade}-macse.log"
    resources:
        mem_mb=config['mem_mb']['macse']
    shell:
        """
        find {input.fasta_dir}/*.{params.fasta_ext} | \
        parallel --jobs {threads} "macse -prog alignSequences \
         -seq {{}} \
         -gc_def {params.gc} \
         -out_NT {output.outdir}/{{/.}}.NT.fasta \
         -out_AA {output.outdir}/{{/.}}.AA.fasta"  
        """

rule concat_align:
    """
    Concatenates alignments
    """

    input:
        aln_dir=directory("output/align/macse/{clade}"),
        aln_done="output/align/macse/{clade}/{clade}.done"
    output:
        aln_cat="output/align/macse/concat/{clade}.aligned.cat.fasta",
        partition="output/align/macse/concat/{clade}.aligned.cat.partitions"
    log:
        "output/logs/align/macse/{clade}-cat.log"
    run:
        aln_files = glob(join(input.aln_dir, '*.NT.fasta'))
        aln_df, len_series = load_alignments(aln_files)
        fasta, partitions, lens = cat_alignments(aln_df, len_series)
        with open(output.aln_cat, 'w') as f:
            f.write(fasta)
        with open(output.partition, 'w') as f:
            f.write(partitions)

