#!/usr/bin/python

import os, math, sys, ROOT, subprocess
ROOT.gROOT.SetBatch(True)

name_list = ["monoWjjIsrDesD1m50", "monoWjjIsrDesD1m1300", "monoWjjIsrDesD5m50", "monoWjjIsrDesD5m1300", "monoWjjIsrDesD9m50", "monoWjjIsrDesD9m1300", "monoWjjIsrConD1m50", "monoWjjIsrConD1m1300", "monoWjjIsrConD5m50", "monoWjjIsrConD5m1300", "monoWjjIsrConD9m50", "monoWjjIsrConD9m1300", "monoWjjWwxxm50", "monoWjjWwxxm1300", "monoWjjSimDesSsdMed100m50", "monoWjjSimDesSsdMed1500m50", "monoWjjSimDesSsdMed1500m1300", "monoWjjSimDesSvdMed100m50", "monoWjjSimDesSvdMed1500m50", "monoWjjSimDesSvdMed1500m1300", "monoWjjSimDesTsdMed100m50", "monoWjjSimDesTsdMed1500m50", "monoWjjSimDesTsdMed1500m1300", "monoWjjSimConSsdMed100m50", "monoWjjSimConSsdMed1500m50", "monoWjjSimConSsdMed1500m1300", "monoWjjSimConSvdMed100m50", "monoWjjSimConSvdMed1500m50", "monoWjjSimConSvdMed1500m1300", "monoZjjIsrD1m50", "monoZjjIsrD1m1300", "monoZjjIsrD5m50", "monoZjjIsrD5m1300", "monoZjjIsrD9m50", "monoZjjIsrD9m1300", "monoZjjZzxxm50", "monoZjjZzxxm1300", "monoZjjSimSsdMed100m50", "monoZjjSimSsdMed1500m50", "monoZjjSimSsdMed1500m1300", "monoZjjSimSvdMed100m50", "monoZjjSimSvdMed1500m50", "monoZjjSimSvdMed1500m1300", "monoZjjSimTsdMed100m50", "monoZjjSimTsdMed1500m50", "monoZjjSimTsdMed1500m1300", "monoHbb_mx1_xdxhDh", "monoHbb_mx65_xdxhDh", "monoHbb_mx1000_xdxhDh", "monoHbb_mx1_xgxFhDh", "monoHbb_mx65_xgxFhDh", "monoHbb_mx1000_xgxFhDh", "monoHbb_mx1_zpzp100", "monoHbb_mx65_zpzp100", "monoHbb_mx1000_zpzp100"]

#name_list = ["monoWjjIsrDesD1m50", "monoWjjWwxxm50", "monoZjjSimSsdMed100m50", "monoWjjSimConSsdMed100m50", "monoZjjSimSsdMed1500m1300"] 

ini_count = 0

for name in name_list:
	cmd = "sed -i 's/place_name/" + name + "/g' FrameworkExe/util/hsg5frameworkReadCxAOD.cxx"
	cmd += "; rc compile"
	cmd += "; wait"
	cmd += "; rm -r ../CxAOD_outputs/" + name 
	cmd += "; hsg5frameworkReadCxAOD ../CxAOD_outputs/" + name + " &> cutflows/cutflow_" + name + ".txt"
	cmd += "; sed -i 's/" + name + "/place_name/g' FrameworkExe/util/hsg5frameworkReadCxAOD.cxx"
	print cmd
#	with open("./autocutflow.txt", "w") as outfile:
#                with open("./autocutflow_err.txt", "w") as outfileerr:
#                        subprocess.call(str(cmd), shell = True, stdout = outfile, stderr = outfileerr)


out_file = open('./yields_5fb_signal.txt', 'w')
for name in name_list:
	got_ID = 0
	results_file =open("./cutflows/cutflow_" + name + ".txt", 'r')
	for line in results_file:
		if 'mc14_13TeV' in line and got_ID == 0:
			j = line.index('mc14_13TeV')
			ID = line[j+11:j+17]
			got_ID = 1
		elif 'mJ' in line:
			i = line.index('mJ')
                        my_yield = line[i+3:]           
                        my_yield_fl = float(my_yield)
                        my_yield_fl *= 2.5
                        my_yield_str = str(my_yield_fl)
		elif 'num_mj' in line:
			raw = line[9:].rstrip()
                        out_file.write(name + " (" + ID + ") : " + raw + " (raw), " + my_yield_str + " (yield)\n")
		else:
			continue
out_file.close()
	



# write
#for i_nrj in nrj :
#  # pick the good table
#  if i_nrj == "8TeV"  : count = count_8
#  if i_nrj == "13TeV" : count = count_13
#  # define file
#  out_file = open('./yields.'+ana+'.'+i_nrj+'.txt', 'w')
#  # write
#  for i_count in count :
#    #print str(i_nrj)+" -> "+i_count[0]+" = "+str(i_count[1])+" "+str(i_count[2])
#    out_file.write(str(i_count[0])+"\t\t"+str(i_count[1])+"\t\t"+str(i_count[2])+"\n")
##  out_file.close()
