# GMSC-mapper

Command line tool to query the Global Microbial smORFs Catalog (GMSC)

## Installation

## Example usage
### Default
1. Input is genome contig sequences.

```bash
python3 main.py -i example.fa -o output
```

2. Input is protein sequences
3. Input is nucleotide sequences

### GMSC database/habitat/taxonomy/quality file path can be assigned on your own
```bash
python3 main.py -i example.fa -o output --db exampledb.dmnd
```
```bash
python3 main.py -i example.fa -o output --habitat ref_habitat.txt 
```
```bash
python3 main.py -i example.fa -o output --taxonomy ref_taxonomy.txt 
```
```bash
python3 main.py -i example.fa -o output --quality ref_quality.txt
```

### Habitat/taxonomy/quality annotation is optional
If you don't want to annotate habitat/taxonomy/quality you can use `--nohabitat`/`--notaxonomy`/`--noquality`.
```bash
python3 main.py -i example.fa -o output --nohabitat
```
```bash
python3 main.py -i example.fa -o output --notaxonomy
```
```bash
python3 main.py -i example.fa -o output --noquality
```
## Example Output
The output folder will contain
- Outputs of smORF prediction (Macrel).
- Complete mapping result table, listing all the hits in GMSC, per smORF.(Default:Diamond/MMseqs)
- Habitat annotation of smORFs.(optional)
- Taxonomy annotation of smORFs.(optional)
- Quality annotation of smORFs.(optional)

## Parameters
* `-i/--input`: path to the input genome contig sequence file (FASTA, possibly .gz/.bz2/.xz compressed).

* `--nt-genes`: path to the input nucleotide sequence file (FASTA, possibly .gz/.bz2/.xz compressed).

* `--aa-genes`: path to the input protein sequence file (FASTA, possibly .gz/.bz2/.xz compressed).

* `-o/--output`: Output directory (will be created if non-existent).

* `--db`: path to the GMSC database file(.dmnd).

* `--habitat`: path to the habitat file.

* `--nohabitat`: use this if no need to annotate habitat.

* `--taxonomy`: path to the taxonomy file.

* `--notaxonomy`: use this if no need to annotate taxonomy.

* `--quality`: path to the quality file.

* `--noquality`: use this if no need to annotate quality.


