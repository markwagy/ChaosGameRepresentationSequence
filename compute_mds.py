from sklearn import manifold
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob

def loadAnnotations():
    all_annotation_fnames = glob.glob('tmp/atn*.pkl')
    # sort annotations by fname to keep them lined up with the order of MDS embedding
    sorted(all_annotation_fnames)
    annotations = [pickle.load(open(fname)) for fname in all_annotation_fnames]
    return annotations

def run(DSSIM_matrix):
    mds = manifold.MDS(n_components=2,
                       max_iter=3000,
                       eps=1e-9,
                       random_state=np.random.RandomState(seed=123),
                       dissimilarity='precomputed', n_jobs=-1)
    print "MDS done"
    results = mds.fit(DSSIM_matrix)
    np.save('MDS_fit', results)
    coords = results.embedding_
    #plt.subplots_adjust(bottom=0.1)
    annotations = loadAnnotations()
    #plt.scatter(coords[:, 0], coords[:, 1], marker='o')
    #plt.show()
    return coords, annotations


if __name__ == '__main__':
    DSSIM_matrix = np.load('DSSIM.npy')
    run(DSSIM_matrix)