import cdo
import multiprocessing as mp
import os


def rotate_and_delete(infile, outfile, debug=True):
    client = cdo.Cdo(debug=debug)
    if client.debug:
        print("Rotating {}...".format(infile.split("/")[-1]))
    client.sellonlatbox("-180,180,-90,90", input=infile, output=outfile)
    if client.debug:
        print("Done. Removing {}...".format(infile.split("/")[-1]))
    os.remove(infile)
    if client.debug:
        print("Done.")


def main():
    months = ["{:02d}".format(m) for m in range(1, 13)]
    years = ["{:02d}".format(y) for y in range(3, 20)]
    yms = [y + m + ".nc" for m in months for y in years]
    yms = ["1302.nc", "1402.nc"]

    args = [
        (
            "/media/zegheim/Justin/nc/cams_gfas_ga_" + ym,
            "/media/zegheim/Justin/nc_rotated/cams_gfas_ga_" + ym,
        )
        for ym in yms
    ]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.starmap_async(rotate_and_delete, args).get()


if __name__ == "__main__":
    main()
