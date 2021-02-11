# AbCluster
"AbCluster" is a tool for antibody sequence analysis including translation, numbering and clustering.

# Requirement
* Ubuntu 16.04
* Python 3.7
* HMMER 3.3
* ANARCI

# Usage
class AbCluster(path)
parameters
* path: a fasta file including antibody DNA or amino acid sequences

class AbCluster.translate(trans, restrict)
parammeters
* trans(True/False): translate the sequences to amino acids when input data are dna sequences
* restrict('ig'/'heavy'/'light'): restrict analysis antibody region

class AbCluster.classification(df, model, n_clusters=8)
parameters
* df: pandas dataframe ordered antibody sequences from "translate" method
* model('kmeans'/'agglomerative_clustering'): algorithm for clustering
* n_clusters[int]: numbers of clustering groups 

# Example

a = AbCluster("~/work/test.fasta")

b = a.translate(trans=True, restrict='ig')

| H1 | H2 | ... | L1 | L2 | ... | germline | name  | species | Hch  | Lch  | Hch+Lch  |
|:--:|:--:|:---:|:--:|:--:|:---:|:--------:|:-----:|:-------:|:----:|:----:|:--------:|
| Q  | V  | ... | D  | I  | ... | IGKV1-27 | A0001 | human   | QV.. | DI.. | QV..DI.. |


c = a.classification(df=b, model='kmeans', n_clusters=8)

| H1 | H2 | ... | L1 | L2 | ... | germline | name  | species | Hch  | Lch  | Hch+Lch  | pred |
|:--:|:--:|:---:|:--:|:--:|:---:|:--------:|:-----:|:-------:|:----:|:----:|:--------:|:----:|
| Q  | V  | ... | D  | I  | ... | IGKV1-27 | A0001 | human   | QV.. | DI.. | QV..DI.. |  1   |
| .  | .  | ... | .  | .  | ... | ........ | A0002 | human   | .... | .... | ........ |  3   |
| .  | .  | ... | .  | .  | ... | ........ | A0003 | human   | .... | .... | ........ |  4   |
| .  | .  | ... | .  | .  | ... | ........ | A0004 | human   | .... | .... | ........ |  6   |
| .  | .  | ... | .  | .  | ... | ........ | A0005 | human   | .... | .... | ........ |  8   |
| .  | .  | ... | .  | .  | ... | ........ | A0006 | human   | .... | .... | ........ |  2   |
| .  | .  | ... | .  | .  | ... | ........ | A0007 | human   | .... | .... | ........ |  3   |
| .  | .  | ... | .  | .  | ... | ........ | A0008 | human   | .... | .... | ........ |  7   |
| .  | .  | ... | .  | .  | ... | ........ | A0009 | human   | .... | .... | ........ |  5   |



# Note

Use Linux for using ANARCI

# Author

* BioSpace
* biospaceinfo.gmail.com

# License

"AbCluster" is under Apache Software License.
