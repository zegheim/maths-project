import cdo
import glob
import multiprocessing as mp


def crop(infile, outfile, extent):
    client = cdo.Cdo()
    extent = ",".join(map(str, extent))
    print("Cropping {} with extent={}".format(infile, extent))
    client.sellonlatbox(extent, input=infile, output=outfile)
    print("Done. Results saved at {}".format(outfile))


def main():
    infiles = glob.glob("/media/zegheim/Justin/nc/*.nc")
    args = [
        (infile, infile.replace("/nc/", "/nc_aus/"), (112.9, 159.2, -55.2, -9.2))
        for infile in infiles
    ]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.starmap_async(crop, args).get()


if __name__ == "__main__":
    main()
