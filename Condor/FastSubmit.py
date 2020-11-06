import os
import time

# fileList = open('../input/fileLabels.txt','r')
# fileIDList = fileList.readlines()
fileIDList = ["QCD17_Pt_600to800\n"]
# for i in range(2):
for i in range(len(fileIDList)): # changing Arguments in submit.jdl
	rfile = open('submit.jdl','r+')
	f1 = rfile.readlines()

	fileID = fileIDList[i]

	f1[5] = "Output = " + fileID[:-1] + ".out\n"
	f1[6] = "Error = " + fileID[:-1] + ".err\n"
	f1[7] = "Log = " + fileID[:-1] + ".log\n"

	f1[12] = "Arguments = " + fileID

	rfile.seek(0)
	rfile.writelines(f1)
	rfile.truncate()
	rfile.close()

	os.system("condor_submit submit.jdl")

	time.sleep(2)
