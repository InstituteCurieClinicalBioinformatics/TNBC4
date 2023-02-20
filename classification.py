from argparse import ArgumentParser
import subprocess
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "scripts"))
from geneExpressionComparison import *

def main(inFile, outDir, mode, dataset, split, threshold, genesListFolder, gitDir):
    code = subprocess.call(f"Rscript {gitDir}/scripts/clustering.R -i {inFile} -g {genesListFolder}/ -o {outDir} -m {mode} -s {split}", shell=True)
    # if code == 0:
        # computeTestStats(os.path.join(outDir, "final_clusters.csv"), dataset, outDir)
        # subprocess.call(f"Rscript {gitDir}/scripts/pathways.R -i {os.path.join(outDir, f'{dataset}_pval.csv')} -o {outDir}/ -t {threshold}", shell=True)

if __name__ == "__main__":
    parser = ArgumentParser(description="-Launch affymetrix analysis")

    parser.add_argument("-d", "--dataset", required=True, help="Dataset ID. Used to name output file.")
    parser.add_argument("-g", "--genesListFolder", required = True, help = "Folder containing gene list for each subtype", default = os.path.join(os.path.dirname(os.path.realpath(__file__)), "genes"))
    parser.add_argument("-i", "--inputFile", required=True, help="Path to the input file.")
    parser.add_argument("-m", "--mode", required=True, help="Data type to cluster. Must be rnaseq or microarray", choices = ["rnaseq", "microarray"])
    parser.add_argument("-o", "--outDir", required=True, help="Path to the output folder.")
    parser.add_argument("-s", "--split", help="Number of data split to perform for crossing gene name with probes name. Depending on specs machine, do not perform split can lead to out of memory exception", type=int, default=200)
    parser.add_argument("-t", "--threshold", help="Pvalue cutoff", type=float, default=0.05)

    args = parser.parse_args()
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    gitDir = os.path.dirname(os.path.realpath(__file__))

    main(args.inputFile, args.outDir, args.mode, args.dataset, args.split, args.threshold, args.genesListFolder, gitDir)
