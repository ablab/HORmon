# HORmon
HORmon is a tool for annotation of alpha satellite arrays in centromeres of a newly assembled human genome. HORmon consists of two modules:
* **Monomer Inference** extract draft human monomers based on the given alpha-satellite consensus template and centromeric sequence.
* **HORmon** polishes monomers extracted on the previous stage to make it consistent with Centromere Evolution postulate, extract HORs and decompose centromeric sequence into HORs.

## Installation
Requirements:
* Python3.6+
  * biopython
  * clustalo
  * joblib
  * python-edlib
  * setuptools
  * networkx
  * pygraphviz
  * stringdecomposer

The required python packages can be installed through conda using
```
conda install --file requirement.txt
```

Installing from source
```
git clone https://github.com/ablab/centromere-architect.git
cd centromere-architect
???
```

Then, HORmon is available as `monomer_inference` and `HORmon` 

## Quick start
### Monomer Inference
Monomer Inference script needs two parameters: (1) (centromeric) sequence and (2) monomer template:

```
python3 src/monomer_inference.py -seq test_data/cen8toy.fa -mon test_data/AlphaSat.fa
```

Resulting monomers can be found in ```final/monomers.fa``` and sequence annotation in ```final/final_decomposition.tsv```.

### HORmon
HORmon takes as input: (1) (centromeric) sequence, (2) monomers from "monomer_inference" stage, (3) centromere id for nicer output, (4) the minimum number of occurrence for monomers, monomers pairs and HOR and (5) output folder
```
python3 src/HORmon.py --seq test_data/cen8toy.fa --mon test_data/cen8_mn.fa --cen-id 8 --monomer-thr 2 --edge-thr 2 --min-traversals 2 -o toy8
```

Output:
* `toy8/mn.fa` -- final monomers
* `toy8/final_decomposition.tsv` -- monomer decomposition
* `toy8/HORs.tsv` -- HORs description
* `toy8/HORdecomposition.tsv` -- HORs decomposition

## Feedback and bug reports
Your comments, bug reports, anf suggestions are very welcomed. They will help up to further improve HORmon. 

You can leave your comments and bug reports at [our GitHub repository tracker](https://github.com/ablab/centromere-architect/issues). 
