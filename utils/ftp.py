import ftputil

with ftputil.FTPHost("ftp.cdc.noaa.gov", "anonymous", "s1545572@sms.ed.ac.uk") as host:
    host.download_if_newer(
        "/Datasets/cpc_global_temp/tmax.2015.nc",
        "/media/zegheim/Justin_SSD/cpc_global_temp/tmax.2015.nc",
    )
    host.download_if_newer(
        "/Datasets/cpc_global_temp/tmin.2015.nc",
        "/media/zegheim/Justin_SSD/cpc_global_temp/tmin.2015.nc",
    )
