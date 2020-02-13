import cdo
import glob

from datetime import datetime


def crop(client, infile, outfile, extent):
    extent = ",".join(map(str, extent))
    client.sellonlatbox(extent, input=infile, output=outfile)


def main():
    client = cdo.Cdo()
    infiles = glob.glob("/media/zegheim/SEAGATE_POR/data_sst_MUR/raw/*.nc")
    start = datetime.now()
    for idx, infile in enumerate(infiles):
        print(
            "Progress: {:0>3d}/{}. Time elapsed: {:0>8.2f}s".format(
                idx, len(infiles), (datetime.now() - start).total_seconds()
            ),
            sep="\r",
            flush=True,
        )
        outfile = infile.replace("/raw/", "/cropped/")
        extent = (-5, -4, 52, 54)
        try:
            crop(client, infile, outfile, extent)
        except:
            print("Could not crop {}. Skipping".format(infile.split("/")[-1]))
            continue

    print(
        "Done. Total time elapsed: {:.2f}s".format(
            (datetime.now() - start).total_seconds()
        )
    )


if __name__ == "__main__":
    main()

