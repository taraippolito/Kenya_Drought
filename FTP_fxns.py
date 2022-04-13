

def download_file(filename):
	import os, sys, os.path
	from ftplib import FTP
    #FTP login
	ftp = FTP("ftp.chc.ucsb.edu")
	ftp.login()
	directory = '/pub/org/chc/products/CHIRTSdaily/v1.0/global_tifs_p05/Tmax/2000/'
	ftp.cwd(directory)

	ddir='F:\\CHIRTS'
	local_filename = os.path.join(ddir, filename)
	print(local_filename)
	file = open(local_filename, 'wb')
	ftp.retrbinary('RETR '+ filename, file.write)
	file.close()
	ftp.quit()