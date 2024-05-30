configfile: "config.yaml"


comps = config["comparisons"].keys()
all_samples = set([s for c in comps for s in config["comparisons"][c]])

print(all_samples)


rule seq_lengths:
    input:
        fas=expand("data/{sample}.fa", sample=all_samples),
    output:
        "results/seq_lengths.csv",
    shell:
        """
        python scripts/seq_lengths.py \
            --fastas {input.fas} \
            --out {output}
        """


rule build_graph:
    input:
        ref=lambda w: f'data/{config["comparisons"][w.comp][0]}.fa',
        qry=lambda w: f'data/{config["comparisons"][w.comp][1]}.fa',
    output:
        pan="results/{comp}/graph.json",
    shell:
        """
        pangraph build \
            --circular \
            -s 20 \
            -a 100 \
            -b 5 \
            {input.ref} {input.qry} \
            > {output.pan}
        """


rule block_stats:
    input:
        pan=rules.build_graph.output,
    output:
        stats="results/{comp}/block_stats.csv",
    shell:
        """
        python scripts/block_stats.py \
            --graph {input.pan} \
            --out {output.stats}
        """


rule export_gfa:
    input:
        pan=rules.build_graph.output,
    output:
        exp=directory("results/{comp}/export"),
    shell:
        """
        pangraph export \
            {input.pan} \
            --no-duplications \
            -ell 0 \
            --output-directory {output.exp}
        """


rule core_alignments:
    input:
        pan=rules.build_graph.output,
    output:
        aln=directory("results/{comp}/core_alignments"),
    shell:
        """
        python scripts/core_blocks_alignments.py \
            --graph {input.pan} \
            --out_fld {output.aln}
        """


rule dotplot:
    input:
        pan=rules.build_graph.output,
        lens="results/seq_lengths.csv",
    output:
        "results/{comp}/dotplot.html",
    shell:
        """
        python scripts/dotplot.py \
            --graph {input.pan} \
            --seq_lengths {input.lens} \
            --output {output}
        """


rule all:
    input:
        expand(rules.block_stats.output, comp=comps),
        expand(rules.export_gfa.output, comp=comps),
        expand(rules.core_alignments.output, comp=comps),
        expand(rules.dotplot.output, comp=comps),
