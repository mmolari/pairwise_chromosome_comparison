# pairwise chromosome comparison

simple pipeline for pairwise chromosome comparison

## setup

### installing the tools

- Install and activate the provided conda environment
```sh
conda env create -n paircomp --file conda_env.yml
```
- install pangraph following [the installation instructions](https://neherlab.github.io/pangraph/#Installation).

### adding the data

- add the chromosomes to compare as fasta files in `data/XXX.fa` and `data/YYY.fa`, where `XXX` and `YYY` stand for the unique file names. Each should contain *a single contig* which should have a simple id without white spaces.
- in the `comparisons` entry of the `config.yaml` file add the desired comparison:
```yaml
comparisons:
  comp_1: ["XXX", "YYY"]
```

## running the pipeline

after activating the conda environment run the pipeline with:
```sh
snakemake -c1 all
```

## results

For any given comparison, the pipeline produces the following results, stored in the `results/{comp_id}` folder.

### block statistics

The `block_stats.csv` file contains a list of all of the blocks in the graph with summary statistics:

|            | count | n. strains |    len | duplicated | core  | category   |
| :--------- | ----: | ---------: | -----: | :--------- | :---- | :--------- |
| OGIXJBTVTV |     2 |          2 | 326946 | False      | True  | core       |
| XEVNWFKICB |     2 |          2 | 224824 | False      | True  | core       |
| IRWDQBJGUJ |     2 |          2 | 224159 | False      | True  | core       |
| YTFFFPIQSD |     1 |          1 |   1418 | False      | False | accessory  |
| UZWIYCEDJP |     1 |          1 |    869 | False      | False | accessory  |
| UJFECAATVC |    23 |          2 |   1055 | True       | False | duplicated |
| NPESBINNJX |    21 |          2 |    142 | True       | False | duplicated |
| NBPXZUNQOA |    14 |          2 |   3243 | True       | False | duplicated |

- the block **count** is the total number of occurrences of the block in the two paths.
- the block **n. of strains** is the number of strains in which the block appears, either one or two.
- the block **length** is the length of the block consensus sequence.
- the **duplicated** flag indicates whether the block occurrs twice in any path.
- the **core** flag indicates block that are found exactly once in each genome.
- the **category** value can be either *core*, *accessory* (for blocks that are found only once in a path but not the other) or *duplicated*


## ToDo

- [x] extract alignments of all core regions
- [x] extract list of mutations: three dataframes: SNPs, insertions and deletions for core blocks.
- [x] dotplot
  - [x] add private segments on the sides (sub-plots)
  - [x] add periodic boundary conditions
- [x] block stats
- [x] core alignments
- [ ] detect non-syntenic duplications by gluing the paths
- [ ] position in the genome for the mutations / indels
- [ ] block list with positions in the two genomes