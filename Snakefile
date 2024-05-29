configfile: "config.yaml"


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


comps = config["comparisons"].keys()


rule all:
    input:
        expand(rules.block_stats.output, comp=comps),
        expand(rules.export_gfa.output, comp=comps),
        expand(rules.core_alignments.output, comp=comps),
