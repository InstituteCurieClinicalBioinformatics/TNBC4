import pandas as pd
from scipy.stats import ranksums
import pickle
from statsmodels.stats.multitest import multipletests
from argparse import ArgumentParser
import os

def pval(clusterisedData, cluster, dataset):
    clusters = pd.read_csv(clusterisedData, delimiter=",", index_col = 0)
    one = clusters[clusters["cluster"] == cluster]
    others = clusters[clusters["cluster"] != cluster]
    pval = []
    for i in range(one.shape[1] - 1):
        pval.append(ranksums(one.iloc[:,i], others.iloc[:,i]).pvalue)
    pickle.dump(pval, open("{}_{}_VS_ALL.pickle".format(dataset, cluster), "wb"))

def pvalCorrection(dataset, cluster):
    pval = pickle.load(open("{}_{}_VS_ALL.pickle".format(dataset, cluster), "rb"))
    correction = multipletests(pval, method="fdr_bh")
    pickle.dump(correction[1], open("{}_{}_VS_ALL_corrected.pickle".format(dataset, cluster), "wb"))

def addCorrection(dataset, clusterisedData, cluster):
    pval = pickle.load(open("{}_{}_VS_ALL_corrected.pickle".format(dataset, cluster), "rb"))
    clusterisedData["{}_VS_ALL".format(cluster)] = pval
    return clusterisedData

def annotClusters(clusterisedData, conclusions, outFile):
    clusterisedData["conclusion"] = clusterisedData["cluster"]
    conclusions = pd.read_csv(conclusions, sep = ",", index_col=0)
    conclusion = conclusion["conclusion"].str.split("=", expand = True)
    conclusion = dict(zip(conclusion[1], conclusion[0]))
    for subtype in conclusion.keys():
        if conclusion[conclusion["conclusion"] == subtype].shape[0] == 1:
            clusterisedData["conclusion"] = clusterisedData["conclusion"].replace(subtype, conclusion[subtype])
    clusterisedData.to_csv(outFile, sep = ",")

def computeTestStats(inFile, dataset, outDir):
    df = pd.read_csv(inFile, delimiter=",", index_col = 0)
    clusters = df.pop("cluster").tolist()
    df = df.T
    for i in range(4):
        pval(inFile, i+1, dataset)
        pvalCorrection(dataset, i+1)
        df = addCorrection(dataset, df, i+1)
        print("{} - cluster {} : done".format(dataset, i+1))
    pickles = [f for f in os.listdir() if "pickle" in f]
    for pickle in pickles:
        os.remove(pickle)
    df.to_csv("{}/{}_pval.csv".format(outDir, dataset), sep=",", index=True)

if __name__ == "__main__":
    parser = ArgumentParser(description="-Script to compute mean difference and Mann-Withney on TNBC data")

    parser.add_argument("-d", "--dataset", required=True, help="Dataset ID. Used to name output file")
    parser.add_argument("-i", "--inpuFile", required=True, help="Path of the file final_clusters.csv")
    parser.add_argument("-o", "--outDir", required=True, help="Path to the output folder")

    args = parser.parse_args()

    computeTestStats(args.inFile, args.dataset, args.outDir)