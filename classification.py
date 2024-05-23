from argparse import ArgumentParser
import subprocess
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "scripts"))
from geneExpressionComparison import *

def main(inFile, outDir, mode, cleaning, nbTopGenes, percentFilter, dataset, split, threshold, genesListFolder, gitDir):
    if mode == "rnaseq" and cleaning:
        subprocess.call(f"python3 {gitDir}/scripts/cleanRNAinput.py -i {inFile} -n {nbTopGenes} -p {percentFilter} -o {outDir}/{dataset}_cleaned.tsv", shell = True)
        subprocess.call(f"Rscript {gitDir}/scripts/clustering.R -i {outDir}/{dataset}_cleaned.tsv -g {genesListFolder}/ -o {outDir} -m {mode} -s {split}", shell=True)
    else:
        subprocess.call(f"Rscript {gitDir}/scripts/clustering.R -i {inFile} -g {genesListFolder}/ -o {outDir} -m {mode} -s {split}", shell=True)
    computeTestStats(os.path.join(outDir, "final_clusters.csv"), dataset, outDir)
    subprocess.call(f"Rscript {gitDir}/scripts/pathways.R -i {os.path.join(outDir, f'{dataset}_pval.csv')} -o {outDir}/ -t {threshold}", shell=True)
    annotClusters(f'{outDir}/final_clusters.csv', f'{outDir}/conclusion.csv', f'{outDir}/final_attribution.csv')
    subprocess.call(f"rm -rf {outDir}/CACHE", shell = True)

if __name__ == "__main__":
    parser = ArgumentParser(description="-Launch affymetrix analysis")

    parser.add_argument("-c", "--cleaning", action='store_true', help="Clean input rnaseq data")
    parser.add_argument("-d", "--dataset", required=True, help="Dataset ID. Used to name output file.")
    parser.add_argument("-g", "--genesListFolder", help = "Folder containing gene list for each subtype", default = os.path.join(os.path.dirname(os.path.realpath(__file__)), "genes"))
    parser.add_argument("-i", "--inputFile", required=True, help="Path to the input file.")
    parser.add_argument("-m", "--mode", required=True, help="Data type to cluster. Must be rnaseq or microarray", choices = ["rnaseq", "microarray"])
    parser.add_argument("-n", "--nbTopGenes", help="Number of genes with the lowest std to keep", default = 1000, type = int)
    parser.add_argument("-o", "--outDir", required=True, help="Path to the output folder.")
    parser.add_argument("-p", "--percentFilter", help="Percent of samples for which gene must have value != 0. Otherwise, it's removed", default = 0.3, type = float)
    parser.add_argument("-s", "--split", help="Number of data split to perform for crossing gene name with probes name. Depending on specs machine, do not perform split can lead to out of memory exception", type=int, default=200)
    parser.add_argument("-t", "--threshold", help="Pvalue cutoff", type=float, default=0.05)

    args = parser.parse_args()
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    gitDir = os.path.dirname(os.path.realpath(__file__))

    main(args.inputFile, args.outDir, args.mode, args.cleaning, args.nbTopGenes, args.percentFilter, args.dataset, args.split, args.threshold, args.genesListFolder, gitDir)
