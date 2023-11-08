#! /usr/bin/python

###########################
#import

import sys, getopt
import os
###########################
#sub functions

def usage():
	print('''This script is used to transform differential results of array to excell.
Usage	python '''+sys.argv[0]+''' [options]
Options:
	-d --dir	directory containing all differential results of each compare.
	-c --comp	compare information file
	-o --output	output directory.
	-v --version	display version information.
	-h --help	display the usage.\n''')
	sys.exit()

#########################
#get options

opts, args = getopt.getopt(sys.argv[1:], "hi:d:c:o:v", ["help", "dir=", "comp=", "output=", "version"])

output = "./"

for op, value in opts:
	if op in ("-d", "--dir"):
		dire = value
	elif op in ("-c", "--comp"):
		comp = os.path.abspath(value)
	elif op in ("-o", "--output"):
		output = value
	elif op in ("-v", "--version"):
		version()
	elif op in ("-h", "--help"):
		usage()

if not ('dire' in dir() and 'comp' in dir()):
	usage()
	sys.exit()

#########################
#main programe
from Yilib import table2excel
from math import log
from glob import glob

def headerbuilt_all(thefiles, header):
	## Get all data for all files
	newdic = {}
	for i in thefiles:
		headertmp = header%(i[:-4])
		newdic[i[:-4]] = [headertmp]
	return newdic

def headerbuilt_DE(thefiles, header, comp):
	## Get statistic data from compare file
	ofile = open(comp, "rU")
	line = ofile.readline().rstrip("\n")
	gpdic = {}
	while line:
		line = ofile.readline().rstrip("\n")
		if line == "":
			break
		sline = line.split("\t")
		gpdic[sline[0] + "_vs_" + sline[1]] = [sline]
	ofile.close()
	## build header for each file
	newdic = {}
	for i in thefiles:
		if i.startswith("up_"):
			logeqtion = ">"
			tagtmp = i[3:-4]
			headertmp = header%(logeqtion, '%.3f'%log(float(gpdic[tagtmp][0][2]), 2), logeqtion, '%.1f'%float(gpdic[tagtmp][0][2]), '%.2f'%float(gpdic[tagtmp][0][3]), '%.2f'%float(gpdic[tagtmp][0][4]), i[:-4])
		elif i.startswith("down_"):
			logeqtion = "<"
			tagtmp = i[5:-4]
			headertmp = header%(logeqtion, '%.3f'%-log(float(gpdic[tagtmp][0][2]), 2), logeqtion, '%.4f'%(1/float(gpdic[tagtmp][0][2])), '%.2f'%float(gpdic[tagtmp][0][3]), '%.2f'%float(gpdic[tagtmp][0][4]), i[:-4])
		elif i.startswith("up+down"):
			logeqtion = ">"
			tagtmp = i[8:-4]
			headertmp = header%(logeqtion, '%.3f'%log(float(gpdic[tagtmp][0][2]), 2)+' or log2FC < '+ '%.3f'%-log(float(gpdic[tagtmp][0][2]), 2), logeqtion, '%.1f'%float(gpdic[tagtmp][0][2]) + ' or FC < '+ '%.4f'%(1/float(gpdic[tagtmp][0][2])), '%.2f'%float(gpdic[tagtmp][0][3]), '%.2f'%float(gpdic[tagtmp][0][4]), i[:-4])
		newdic[i[:-4]] = [headertmp]
	return newdic

compo=open(comp, 'rU')
cc=compo.readlines()
compo.close()
all_diff=[]
diff_file=[]
for i in range(1,len(cc)):
	lines=cc[i]
	lines=lines.strip('\n').split('\t')
	all_di=lines[0]+'_vs_'+lines[1]+'.txt'
	all_diff.append(all_di)
	diff_up='up_'+lines[0]+'_vs_'+lines[1]+'.txt'
	diff_down='down_'+lines[0]+'_vs_'+lines[1]+'.txt'
	updown='up+down_'+lines[0]+'_vs_'+lines[1]+'.txt'
	diff_file.append(diff_up)
	diff_file.append(diff_down)
	diff_file.append(updown)

os.chdir(dire)
thefiles = all_diff
header='''#Condition pair: %s'''
header = headerbuilt_all(thefiles, header)
table2excel(thefiles, "All_Differential_result.xlsx", header)

os.chdir(dire+'/sep')
thefiles = diff_file
header = '''#log2FC %s %s
#FC %s %s
#p-value cut-off <= %s
#q-value cut-off <= %s
#Condition pair: %s'''
header = headerbuilt_DE(thefiles, header, comp)
table2excel(thefiles, "Differential_result.xlsx", header)




