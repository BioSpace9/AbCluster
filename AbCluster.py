import os
import shutil
import subprocess
import re
from collections import OrderedDict
import copy

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from sklearn import metrics
from sklearn.cluster import KMeans

class AbCluster:
    
    def __init__(self, path):
        
        self.path = path
        self.result_path = '{}/result'.format(os.path.dirname(self.path))
        os.mkdir(self.result_path)
    
    def translate(self, trans, restrict):

        if trans == True:
            self._dna_to_aa()
        else:
            shutil.copy(self.path, "{}/trans.fasta".format(self.result_path))
            
        df = self._anarci(restrict)

        return df
    
    
    def classification(self, df, model, n_clusters=8):
        
        if model == 'kmeans':
            self._kmeans(df, n_clusters=n_clusters)
            
        elif model == 'agglomerative_clustering':
            self._agglomerative_clustering(df, n_clusters=n_clusters)
    
        
    def _orfs_trans(self, seq, trans_table, min_aa_length):
        
        out = []
        seq_len = len(seq)
        for nuc in [seq]:
            for frame in range(3):
                trans = str(nuc[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end - aa_start >= min_aa_length:
                        start = frame + aa_start * 3 + 1
                        end = min(seq_len, frame + aa_end * 3 + 3)
                        length = int((end - start) / 3 - 1)
                        out.append((length, start, end, trans[aa_start:aa_end]))
                    aa_start = aa_end + 1
        out.sort()
        return out
        
    
    def _dna_to_aa(self):
        
        seq_records = list(SeqIO.parse(self.path, "fasta"))
        table = 11
        min_aa_len = 10
        
        aa_records = seq_records
        for aa_record in aa_records:
            aa_tmp = self._orfs_trans(aa_record.seq, table, min_aa_len)
            if not aa_tmp:
                aa_record.seq = Seq("")
            else:
                aa_record.seq = Seq(aa_tmp[-1][3])
                
        SeqIO.write(aa_records, '{}/trans.fasta'.format(self.result_path), 'fasta')
    
    
    def _sortedStringList(self, array=[]):
        
        sortDict = OrderedDict()
        for splitList in array:
            sortDict.update({splitList:[int(x) for x in re.split("(\d+)",splitList)if bool(re.match("\d*",x).group())]})

        return [sortObjKey for sortObjKey,sortObjValue in sorted(sortDict.items(), key=lambda x:x[1])]
    
    
    def _sortedHeadList(self, array=[]):
        
        sortDict=OrderedDict()
        for splitList in array:
            sortDict.update({splitList:[splitList[0]]})

        return [sortObjKey for sortObjKey,sortObjValue in sorted(sortDict.items(), key=lambda x:x[1])]
    
    def _get_colindex(self, df):
    
        h_start = -1
        h_count = 0
        l_start = -1
        l_count = 0

        for idx, column in enumerate(list(df.columns)):
            if column[0] == 'H':
                h_count += 1

                if h_start == -1:
                    h_start = idx

            elif column[0] == 'L':
                l_count += 1

                if l_start == -1:
                    l_start = idx

        h_end = h_start + h_count - 1
        l_end = l_start + l_count - 1

        return h_start, h_end, l_start, l_end    

    def _anarci(self, restrict):
        res = subprocess.call('ANARCI -i /result/trans.fasta -o numbered_sequences.anarci -s kabat -r {} --assign_germline'
                              .format(restrict))  
        
        f = open('{}/numbered_sequences.anarci'.format(self.result_path), "r", encoding="utf-8")
        anarci = f.readlines()

        records = []
        record = {}
        idx = ''
        for line in anarci:
            if idx == '':
                record['name'] = line[1:].split('|')[0]
                idx = 'after_name'

            elif line[:17] == '#|species|v_gene|':
                idx = 'before_germline'

            elif idx == 'before_germline':
                record['species'] = line.split('|')[1]
                record['germline'] = line.split('|')[2]
                idx = 'after_germline'

            elif line[0] == 'H' and (restrict == 'ig' or restrict == 'heavy'):
                record[line[:9].replace(' ', '')] = line[10]
            
            elif line[0] == 'L' and (restrict == 'ig' or restrict == 'light'):
                record[line[:9].replace(' ', '')] = line[10]            

            elif line == '//\n':
                records.append(record)
                idx = ''
                record = {}

        df = pd.DataFrame(records)
        df = df.dropna(subset=['species'])
        df = df.fillna('-')

        recolumns = self._sortedStringList(list(df.columns))
        recolumns = self._sortedHeadList(recolumns)
        df = df.reindex(columns = recolumns)

        df_index = self._get_colindex(df)

        if restrict == 'ig':
            df['Hch'] = ''
            df['Hch+Lch'] = ''
            
        elif restrict == 'heavy':
            df['Hch'] = ''
            
        elif restrict == 'light':
            df['Lch'] = ''

            for i in range(len(df)):
                if (restrict == 'ig') or (restrict == 'heavy'):
                    h_seq = ''
                    for j in range(df_index[0], df_index[1] + 1):
                        h_seq = h_seq + df.iloc[i, j]
                    df['Hch'].iloc[i] = h_seq    
                
                if (restrict == 'ig') or (restrict == 'light'):
                    l_seq = ''
                    for j in range(df_index[2], df_index[3] + 1):
                        l_seq = l_seq + df.iloc[i, j]
                    df['Lch'].iloc[i] = l_seq
                
                if (restrict == 'ig'):
                    df['Hch+Lch'].iloc[i] = h_seq + l_seq
        
        df.to_csv('{}/numbered_sequence.csv'.format(self.result_path), index=False)
        
        return df
        
    def _separate_columns(self, df):
        
        dflist = list(df.columns)
        Xlist = copy.copy(dflist)
        Xlist.remove('species')
        Xlist.remove('germline')
        
        df = pd.get_dummies(df, columns=Xlist)
        
        X = df[Xlist]
        y = df[['species', 'germline']]
        
        return X, y
        
    
    def _kmeans(self, df, n_clusters):
        
        X, y =_separate_columns(df)

        n_init = 10
        tol = 0.0001
        max_iter = 300
        random_state = 1
        
        kmeans = KMeans(n_clusters=n_clusters, n_init=n_init, tol=tol, max_iter=max_iter, random_state=random_state)
        kmeans.fit(X)
        pred = kmeans.predict(X)
        
        df['pred'] = pred
        
        df.to_csv('{}/result_kmeans.csv'.format(self.result_path), index=False)
        
        return df
    
    def _aggclust(self, df, n_clusters):
        
        X, y =_separate_columns(df)

        affinity = 'euclidean'
        linkage = 'ward'
        
        cluster = AgglomerativeClustering(n_clusters =n_clusters, affinity=affinity, linkage=linkage)
        pred = cluster.fit_predict(X)
        
        df['pred'] = pred
        
        df.to_csv('{}/result_aggclust.csv'.format(self.result_path, index=False))

