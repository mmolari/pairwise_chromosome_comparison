configfile: "config.yaml"


comps = config["comparisons"].keys()


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


rule seq_lengths:
    input:
        rules.build_graph.input,
    output:
        "results/{comp}/seq_lengths.csv",
    shell:
        """
        python scripts/seq_lengths.py \
            --fastas {input} \
            --out {output}
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
        lengths=rules.seq_lengths.output,
    output:
        "results/{comp}/dotplot.html",
    shell:
        """
        python scripts/dotplot.py \
            --graph {input.pan} \
            --seq_lengths {input.lengths} \
            --output {output}
        """


rule block_positions:
    input:
        pan=rules.build_graph.output,
    output:
        "results/{comp}/block_positions.csv",
    shell:
        """
        python scripts/block_positions.py \
            --graph {input.pan} \
            --out {output}
        """


rule mutations_positions:
    input:
        pan=rules.build_graph.output,
        lengths=rules.seq_lengths.output,
        alns=rules.core_alignments.output,
        bpos=rules.block_positions.output,
    output:
        "results/{comp}/mutations_positions.csv",
    shell:
        """
        python scripts/mutations_positions.py \
            --graph {input.pan} \
            --lengths {input.lengths} \
            --core_alignments {input.alns} \
            --block_positions {input.bpos} \
            --out_csv {output}
        """


rule minimal_synteny_units:
    input:
        pan=rules.build_graph.output,
    output:
        "results/{comp}/msu/minimal_synteny_units.csv",
    shell:
        """
        python scripts/synteny_units.py \
            --graph {input.pan} \
            --out {output}
        """


rule all:
    input:
        expand(rules.block_stats.output, comp=comps),
        expand(rules.export_gfa.output, comp=comps),
        # expand(rules.core_alignments.output, comp=comps),
        expand(rules.dotplot.output, comp=comps),
        # expand(rules.block_positions.output, comp=comps),
        expand(rules.mutations_positions.output, comp=comps),
        expand(rules.minimal_synteny_units.output, comp=comps),
