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

c = a.classification(df=b, model='kmeans', n_clusters=8)


# Note

Use Linux for using ANARCI

# Author

* BioSpace
* biospaceinfo.gmail.com

# License

"AbCluster" is under Apache Software License.
