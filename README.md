# HORmon

[![Anaconda-Server Badge](https://anaconda.org/bioconda/hormon/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

HORmon is a tool for annotation of alpha satellite arrays in centromeres of a newly assembled human genome. HORmon consists of two modules:
* **Monomer Inference** extracts draft human monomers based on the given alpha-satellite consensus template and centromeric sequence.
* **HORmon** polishes monomers extracted on the previous stage to make it consistent with Centromere Evolution postulate, extract HORs and decompose centromeric sequence into HORs.

HORmon has been used to infer monomers from the recently announced complete human genome assembly of the CHM13 cell line generated by the [Telomere-to-Telomere Consortium](https://github.com/marbl/CHM13#v11).
The data generated in the paper that describes HORmon (Kunyavskaya et al., 2021) can be found at [Figshare](https://figshare.com/articles/dataset/HORmon/16755097/1).
The data includes extracted monomers and HORs from all live alpha satellite arrays in CHM13 cell line, as well as annotations of these arrays.
Jupyter notebook that reproduces figures of the HOR paper is available at [github](https://github.com/TanyaDvorkina/hormon_paper/blob/dev/HORmon_paper.ipynb).

## Installation
Requirements:
* Linux only. Mac OS is not yet supported
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
conda install --file requirements.txt
```

Installing from source
```
git clone https://github.com/ablab/centromere-architect.git
cd centromere-architect
python3 setup.py install --record hormon_files.txt
```

Then, HORmon is available as `monomer_inference` and `HORmon`

Afterward, to uninstall HORmon please run
```
xargs rm -rf < hormon_files.txt
```

## Quick start
### Monomer Inference
Monomer Inference script needs two parameters: (1) (centromeric) sequence and (2) monomer template:

```
monomer_inference -seq test_data/cen8toy.fa -mon test_data/AlphaSat.fa -o toy8_mi
```

Resulting monomers can be found in ```toy8_mi/final/monomers.fa``` and sequence annotation in ```toy8_mi/final/final_decomposition.tsv```.

### HORmon
HORmon takes as input: (1) (centromeric) sequence, (2) monomers from "monomer_inference" stage, (3) centromere id for nicer output, (4) the minimum number of occurrence for monomers, monomers pairs and HOR and (5) output folder
```
HORmon --seq test_data/cen8toy.fa --mon toy8_mi/final/monomers.fa --cen-id 8 --monomer-thr 2 --edge-thr 2 --min-traversals 2 -o toy8
```

Output:
* `toy8/mn.fa` -- final monomers
* `toy8/final_decomposition.tsv` -- monomer decomposition
* `toy8/HORs.tsv` -- HORs description
* `toy8/HORdecomposition.tsv` -- HORs decomposition


***Note***: *HORmon should be launch on each centromere independently. We currently cannot guarantee adequate results, in the case of running HORmon on all centromeres simultaneously.*

## CentromereArchitect 

CentromereArchitect (early version of HORmon), as it is described in the paper, is available at the branch [centromere-architect](https://github.com/ablab/centromere-architect/tree/centromere-architect).
Please cite [Dvorkina et al., 2021](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i196/6319687).

## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcomed.  
Please leave them at [our GitHub repository tracker](https://github.com/ablab/centromere-architect/issues).

# Cite

- [Dvorkina, T., Kunyavskaya, O., Bzikadze, A. V., Alexandrov, I., & Pevzner, P. A. (2021). CentromereArchitect: inference and analysis of the architecture of centromeres. Bioinformatics, 37(Supplement_1), i196-i204.](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i196/6319687)
- Kunyavskaya, O., Dvorkina, T, Bzikadze, A. V., Alexandrov, I., & Pevzner, P. A. (2021). HORmon: automated annotation of human centromeres. bioRxiv
