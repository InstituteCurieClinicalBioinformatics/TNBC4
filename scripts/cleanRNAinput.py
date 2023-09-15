import pandas as pd
from argparse import ArgumentParser

def removeGenes(row):
	if row.tolist().count(0) / len(row.tolist()) > 0.3:
		return False
	return True

def filterOnStandardError(df, outFile):
	df["Standard_error"] = df.std(axis = 1) / len(df.columns.tolist())
	df = df.sort_values(by = ["Standard_error"], ascending = False)
	df.iloc[:1000].to_csv(outFile, sep = "\t")

def main(rnaSeq, outFile):
	df  = pd.read_csv(rnaSeq, sep = "\t", index_col = 0)
	df["Keep"] = df.apply(lambda row: removeGenes(row), axis = 1)
	df = df[df["Keep"] == True]
	filterOnStandardError(df, outFile)

if __name__ == "__main__":
	parser = ArgumentParser(description="Clean RNASeq data")

	parser.add_argument("-i", "--inFile", required=True, help="Path of the RNASeq data")
	parser.add_argument("-o", "--outFile", required=True, help="Path to the clean RNASeq data")

	args = parser.parse_args()

	main(args.inFile, args.outFile)
