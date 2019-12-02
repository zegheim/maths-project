import datetime
import os
from ecmwfapi import ECMWFDataServer

PARAMS = ["80.210", "82.210", "87.210", "88.210", "99.210"]


def get_last_date(d):
    """
    Given a date, returns the date corresponding to the earlier of the last day of the same month/year as the given date, or yesterday.
    """

    last_date = datetime.date(
        year=d.year + d.month // 12, month=d.month % 12 + 1, day=1
    ) - datetime.timedelta(days=1)

    return min(datetime.date.today() - datetime.timedelta(days=1), last_date)


def get_date_ranges(start):
    """
    Provides a list of monthly start-end date tuples for further use.
    """
    end = get_last_date(start)
    date_ranges = [(start, end)]

    while end < datetime.date.today() - datetime.timedelta(days=1):
        start = end + datetime.timedelta(1)
        end = get_last_date(start)
        date_ranges.append((start, end))

    return date_ranges


def retrieve(server, date_range):
    print(
        "Downloading data from {} to {}".format(
            date_range[0].strftime("%Y-%m-%d"), date_range[1].strftime("%Y-%m-%d")
        )
    )

    server.retrieve(
        {
            "class": "mc",
            "dataset": "cams_gfas",
            "date": "/to/".join(d.strftime("%Y-%m-%d") for d in date_range),
            "expver": "0001",
            "levtype": "sfc",
            "param": "/".join(PARAMS),
            "step": "0-24",
            "stream": "gfas",
            "target": "/media/zegheim/Justin/nc/cams_gfas_ga_{}.nc".format(
                date_range[0].strftime("%y%m")
            ),
            "format": "netcdf",
            "time": "00",
            "type": "ga",
            "use": "infrequent",
        }
    )


def main():
    server = ECMWFDataServer()
    start = datetime.date(2003, 1, 1)
    date_ranges = get_date_ranges(start)

    for date_range in date_ranges:
        retrieve(server, date_range)


if __name__ == "__main__":
    main()
