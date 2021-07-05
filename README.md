# CentromereArchitect

## Version 0.5

CentromereArchitect (CA) is the first tool for annotation of alpha satellite arrays in centromeres of a newly assembled human genome.
CA consists of two modules: 
- Monomer Inference allows extraction of human monomers based on the given alpha-satellite consensus template and centromeric sequence.
- HOR Inference allows extraction of HORs from the centromeric sequences using the inferred monomers.


## Installation

Requirements:
- Python3.5+
    - [biopython](https://biopython.org/wiki/Download)
    - [clustalo](https://pypi.org/project/clustalo/)
    - [python-edlib](https://pypi.org/project/edlib/)
    - [joblib](https://pypi.org/project/joblib/)
    - [setuptools](https://pypi.org/project/setuptools/)
    - [StringDecomposer](https://github.com/ablab/stringdecomposer)

## Quick start

### Monomer Inference

Monomer Inference script needs two 1) parameters (centromeric) sequence and 2) monomer template:

```
python3 src/monomer_inference.py -seq test_data/cenXtoy.fasta -mon test_data/AlphaSat.fa
```

Resulting monomers can be found in ```final/monomers.fa``` and sequence annotation in ```final/final_decomposition.tsv```.


### HOR Inference

HOR Inference script needs four parameters 1) (centromeric) sequence, 2) monomers, 3) sequence annotation, and 4) output file name:

```
python src/extract_hors.py test_data/cenXtoy.fasta final/monomers.fa final/final_decomposition.tsv final/hor_decomposition.tsv
```

Resultsing HOR annotation can be found in ```final/hor_decomposition.tsv```.


## Contact

In case of any issues please email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)