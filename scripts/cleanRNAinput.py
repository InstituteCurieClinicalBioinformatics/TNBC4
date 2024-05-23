import pandas as pd
from argparse import ArgumentParser

def removeGenes(row, threshold):
    if row.tolist().count(0) / len(row.tolist()) > threshold:
        return False
    return True

def filterOnStandardError(df, topNb, outFile):
    df["Standard_error"] = df.std(axis = 1) / len(df.columns.tolist())
    df = df.sort_values(by = ["Standard_error"], ascending = False)
    df = df.drop(["Keep", "Standard_error"], axis = 1)
    df.iloc[:topNb].to_csv(outFile, sep = "\t")

def main(rnaSeq, percentFilter, nbTopGenes, outFile):
    df  = pd.read_csv(rnaSeq, sep = "\t", index_col = 0)
    df["Keep"] = df.apply(lambda row: removeGenes(row, percentFilter), axis = 1)
    df = df[df["Keep"] == True]
    filterOnStandardError(df, nbTopGenes, outFile)

if __name__ == "__main__":
    parser = ArgumentParser(description="Clean RNASeq data")

    parser.add_argument("-i", "--inFile", required=True, help="Path of the RNASeq data")
    parser.add_argument("-o", "--outFile", required=True, help="Path to the clean RNASeq data")
    parser.add_argument("-n", "--nbTopGenes", help="Number of genes with the lowest std to keep", default = 1000, type = int)
    parser.add_argument("-p", "--percentFilter", help="Percent of samples for which gene must have value != 0. Otherwise, it's removed", default = 0.3, type = float)

    args = parser.parse_args()

    main(args.inFile, args.percentFilter, args.nbTopGenes, args.outFile)
