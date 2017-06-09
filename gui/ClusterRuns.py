import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

#file = "obs-index_t1_mod.txt"

#run, ra_read, dec_read, glon_read, glat_read = np.loadtxt(file, usecols=(0,1,2,3,4), unpack='True')

def MakeCluster(run,ra_read,dec_read):
  centers  = []

  for ii in range(len(run)):
        gg = 0
        if ra_read[ii] > 180:
            gg= ra_read[ii]-360
        else:
            gg= ra_read[ii]
        #print(glon_read[ii]," -> ",gg)
        centers.append([gg,dec_read[ii]])
   
         
# n_samples put to correspond to the number of runs in the file
# shuffle = 0, otherwise it randomize the positions...
  X, labels_true = make_blobs(n_samples=len(run), centers=centers, cluster_std=0.001, shuffle=0,
                            random_state=0)

  print('Length of X array = %i' % len(X))

# Compute DBSCAN
# eps = 2.0 should be the search radius
  db = DBSCAN(eps=2.0, min_samples=1, metric='euclidean').fit(X)
  core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
  core_samples_mask[db.core_sample_indices_] = True
  labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
  n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

  print('Estimated number of clusters: %d' % n_clusters_)
  print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
  print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
  print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))

##############################################################################
# Print results and make run lists
  from collections import defaultdict
  print_cluster = 42

  runlists = defaultdict(list)

  for j in range(len(run)):
    #if glon_read[j] < 188: 
    #    if glon_read[j] > 187:
    #        print("in loop => ", labels[j], " ", run[j], " ", glon_read[j], " ", glat_read[j]) 
    runlists[labels[j]].append(int(run[j]))

  return runlists
