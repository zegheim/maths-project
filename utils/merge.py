import cdo
import glob


def merge(infile1, infile2, outfile):
    client = cdo.Cdo()
    print("Merging {} and {}".format(infile1, infile2))

