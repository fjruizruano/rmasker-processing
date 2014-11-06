#!/usr/bin/python
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from subprocess import call
import sys
import operator

#usage
#read_rm.py file.fas file_db.fas 10

try:
	data_file = sys.argv[1]
except:
	data_file = raw_input("Introduce FASTA: ")

try:
	file = data_file + ".out"
	open(file)
	
except:
	print "Running RepeatMasker"
	try:
		library = sys.argv[2]
		thread = sys.argv[3]
	except:
		library = raw_input("Introduce DB FASTA file: ")
		thread = raw_input("Introduce number of threads: ")
	call("RepeatMasker -par %s -nolow -no_is -engine crossmatch -lib %s %s" % (thread, library, data_file), shell=True)

thr = 11

data = SeqIO.parse(open(data_file) ,"fasta")
file = open(file).readlines()

######################################################################
test = open("test.fas","w")
test_fun = open("test_fun.fas","w")
def check_orf(secu, begin_d, left_d,sense):
	secu = Seq(secu)
	secu_inv = secu.reverse_complement()
	f1 = secu
	f2 = secu[1:]
	f3 = secu[2:]
	f4 = secu_inv
	f5 = secu_inv[1:]
	f6 = secu_inv[2:]

	frame = ""
	if sense == "+":
		fr = begin_d%3
		if fr == 1:
			frame = f1
		elif fr == 2:
			frame = f3
		elif fr == 0:
			frame = f2

	elif sense == "C":
		fr = left_d%3
		if fr == 1:
			frame = f4
		elif fr == 2:
			frame = f6
		elif fr == 0:
			frame = f5
	test.write(">%s\n%s\n" % (str(begin_d)+"_"+str(left_d)+sense,frame))
	frame_nuc = frame
	frame = frame.translate()
	check = False
	if frame[0] != "*" and frame[-1] != "*":
		check = "*" not in frame
	elif frame[0] == "*":
		check = "*" not in frame[1:]
	elif frame[-1] == "*":
		check = "*" not in frame[0:-1]

	if check == True:
		test_fun.write(">%s\n%s\n" % (str(begin_d)+"_"+str(left_d)+sense,frame_nuc))

	return check
	
######################################################################

print "Creating database of sequences"
dict_seq = {}
for secu in data:
	dict_seq[str(secu.id)] = str(secu.seq)
print "DONE"

# dictionaries
dict_count = {} # counter of sequences
dict_diver = {} # divergences
dict_allen = {} # accumulator of weight

dict_number = {"FULL":0, "EXT3":0, "EXT5":0, "DEL":0, "INS":0, "DEL_INV":0, "INS_INV":0} # count orf sequences
dict_pond = {"FULL":0, "EXT3":0, "EXT5":0, "DEL":0, "INS":0, "DEL_INV":0, "INS_INV":0} # count orf sequences
dict_trans = {"FULL":0, "EXT3":0, "EXT5":0} # count orf sequences
dict_trans_nuc = {"FULL":0, "EXT3":0, "EXT5":0} # count orf sequences
number_trunc_5p = 0
number_trunc_3p = 0
number_trunc_5p_fun = 0
number_trunc_3p_fun = 0
number_trunc_5pp = 0
number_trunc_3pp = 0

match_reads = {} # counter of alignments

# count number of alignment in each read
for line in file[3:]:
	c = -1
	text = line.split()
	for elem in text:
		c += 1
		if elem[0] == "(":
			elem = elem.split("(",1)[1].split(")")[0]
			text[c] = elem
			
	if text[-1] != "*":
#	if text[-1] != "***":
		look = text[4] in match_reads
		if look == True:
			match_reads[text[4]].append(text)
		elif look == False:
			match_reads[text[4]] = [text]

out = open(data_file + ".ful", "w")
outnf = open(data_file + ".nof", "w")

extremes = open(data_file + ".ext.fas", "w")

fulls = open(data_file + ".full.fas", "w")
fulls_fun = open(data_file + ".full.fun.fas", "w")
fulls_def = open(data_file + ".full.def.fas", "w")

insertions = open(data_file + ".insertions.fas", "w")

truncations = open(data_file + ".trunc_nonrte.fas", "w")

eli_def = 0
eli_fun = 0
nuc = 0

#print match_reads

for el in match_reads:
	data = match_reads[el]
	if len(data) == 1:

		#count functional
		data = data[0]

		name = data[4]

		begin_r = int(data[5])
		end_r = int(data[6])+1
		left_r = int(data[7])

		sense = data[8]

		begin_d = int(data[11])
		end_d = int(data[12])
		left_d = int(data[13])

		secu = dict_seq[name]

		nuc += end_r-begin_r #contador nucleotidos
		if sense == "+":
			if begin_d + left_d < 2:
				sequen = secu[begin_r-1:end_r-1]
				fulls.write(">FOR_%s\n%s\n" % (name,sequen))

				check = check_orf(sequen, 1, 1, sense)
				if check == True and len(sequen) == 178:
					eli_fun += 1
					fulls_fun.write(">FOR_%s\n%s\n" % (name,sequen))
				else:
					eli_def += 1
					fulls_def.write(">FOR_%s\n%s\n" % (name,sequen))

		elif sense == "C":
			if begin_d + left_d < 2:
				sequen = secu[begin_r-1:end_r-1]
				sequen = Seq(sequen)
				sequen = str(sequen.reverse_complement())
				fulls.write(">REV_%s\n%s\n" % (name,sequen))

				check = check_orf(sequen, 1, 1, "+")
				if check == True and len(sequen) == 178:
					eli_fun += 1
					fulls_fun.write(">REV_%s\n%s\n" % (name,sequen))
				else:
					eli_def += 1
					fulls_def.write(">REV_%s\n%s\n" % (name,sequen))

		#start classification
		if begin_r + left_r < thr:
			match_reads[el][0] = ["FULL"]+data+["DEF"]
			dict_number["FULL"] += 1
			check = dict_seq[name]
			check = check_orf(check[begin_r-1:end_r-1], begin_d, left_d, sense)
			if check == True:
				dict_trans["FULL"] += 1
				dict_trans_nuc["FULL"] += end_r-begin_r
				match_reads[el][0] = ["FULL"]+data+["FUN"]
			data = "\t".join(match_reads[el][0])
			out.write(data + "\n")

		else:
			thr_sp = 0
			if sense == "+":

				if begin_d < thr:
					match_reads[el][0] = ["EXT5"]+data+["DEF"]
					dict_number["EXT5"] += 1

					secu = dict_seq[name]
#					check = check_orf(secu[begin_r:end_r])
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True:
						dict_trans["EXT5"] += 1
						dict_trans_nuc["EXT5"] += end_r-begin_r
						match_reads[el][0] = ["EXT5"]+data+["FUN"]
					extremes.write(">%s\n%s\n" % ("EXT5p"+str(dict_number["EXT5"]),secu[0:begin_r] ))
					data = "\t".join(match_reads[el][0])
					out.write(data + "\n")

				elif left_d <= thr:
					match_reads[el][0] = ["EXT3"]+data+["DEF"]
					dict_number["EXT3"] += 1
					secu = dict_seq[name]
#					check = check_orf(secu[begin_r:end_r])
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True:
						dict_trans["EXT3"] += 1
						dict_trans_nuc["EXT3"] += end_r-begin_r
						match_reads[el][0] = ["EXT3"]+data+["FUN"]
					extremes.write(">%s\n%s\n" % ("EXT3p"+str(dict_number["EXT3"]), secu[end_r:] ))
					data = "\t".join(match_reads[el][0])
					out.write(data + "\n")

				else:
					match_reads[el][0] = ["DEL"]+data+["DEF"]
					data = "\t".join(match_reads[el][0])
					outnf.write(data + "\n")
					dict_number["DEL"] += 1

				trunc_seq = "NOPE"
				if left_r < thr and begin_r > thr:
					number_trunc_5p += 1
					if left_d < thr_sp:
						number_trunc_5pp += 1
					trunc_seq = secu[0:begin_r]
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >= 0:
						number_trunc_5p_fun += 1
				elif left_r >= thr and left_d < thr and begin_r > thr:
					number_trunc_5p += 1
					trunc_seq = secu[0:begin_r]
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >= 0:
						number_trunc_5p_fun += 1

				if begin_r < thr and left_r > thr:
					number_trunc_3p += 1
					if begin_d < thr_sp:
						number_trunc_3pp += 1
					trunc_seq = secu[end_r:]
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >= 0:
						number_trunc_3p_fun += 1
				elif begin_r >= thr and left_r > thr and begin_d > thr:
					number_trunc_3p += 1
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >= 0:
						number_trunc_3p_fun += 1

				if trunc_seq != "NOPE":
					truncations.write(">%s\n%s\n" % (name, trunc_seq))
					 
			elif sense == "C":

				if begin_d < thr:  # 20:
					match_reads[el][0] = ["EXT3"]+data+["DEF"]
					dict_number["EXT3"] += 1
					secu = dict_seq[name]
#					check = check_orf(secu[begin_r:end_r])
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True:
						dict_trans["EXT3"] += 1
						dict_trans_nuc["EXT3"] += end_r-begin_r
						match_reads[el][0] = ["EXT3"]+data+["FUN"]
					extremes.write(">%s\n%s\n" % ("EXT3m"+str(dict_number["EXT3"]), secu[:begin_r] ))
					data = "\t".join(match_reads[el][0])
					out.write(data + "\n")
				
				elif left_d <= thr: # 20:
					match_reads[el][0] = ["EXT5"]+data+["DEF"]
					dict_number["EXT5"] += 1
					secu = dict_seq[name]
#					check = check_orf(secu[begin_r:end_r])
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True:
						dict_trans["EXT5"] += 1
						dict_trans_nuc["EXT5"] += end_r-begin_r
						match_reads[el][0] = ["EXT5"]+data+["FUN"]
					extremes.write(">%s\n%s\n" % ("EXT5m"+str(dict_number["EXT5"]),secu[end_r:] ))
					data = "\t".join(match_reads[el][0])
					out.write(data + "\n")
				else:
					match_reads[el][0] = ["DEL"]+data+["DEF"]
					data = "\t".join(match_reads[el][0])
					outnf.write(data + "\n")
					dict_number["DEL"] += 1

				trunc_seq = "NOPE"
				if begin_r < thr and left_r > thr:
					number_trunc_5p += 1
					if begin_d < thr_sp:
						number_trunc_5pp += 1
					trun_seq = secu[end_r:]
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >=0:
						number_trunc_5p_fun += 1
				elif begin_r >= thr and left_r > thr and begin_d > thr:
					number_trunc_5p += 1
					trun_seq = secu[end_r:]
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >=0:
						number_trunc_5p_fun += 1
				
				if left_r < thr and begin_r > thr:
					number_trunc_3p += 1
					if left_d < thr_sp:
						number_trunc_3pp += 1
					trunc_seq = secu[0:begin_r]	
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >=0:
						number_trunc_3p_fun += 1
				elif left_r >= thr and left_d < thr and begin_r > thr:
					number_trunc_3p += 1
					trunc_seq = secu[0:begin_r]	
					check = check_orf(secu[begin_r-1:end_r-1], begin_d, left_d, sense)
					if check == True and len(secu[begin_r-1:end_r-1]) >=0:
						number_trunc_3p_fun += 1

				if trunc_seq != "NOPE":
					truncations.write(">%s\n%s\n" % (name, trunc_seq))

	# if there is more than two alignments by sequence	
	else:
		name = data[0][4]
		secu = dict_seq[name]
		c = -1

		classi = ""

		dataone = data[0]
		datatwo = data[1]

		begin_r = int(datatwo[5])
		end_r = int(dataone[6])

		if data[0][8] == data[1][8]:
			if begin_r - end_r < thr:
				classi = "DEL"
				dict_number[classi] += 1
			elif begin_r - end_r >= thr:
				classi = "INS"
				dict_number[classi] += 1
				insertions.write(">INS_%s\n%s\n" % (name,secu[end_r:begin_r]))
			
		elif data[0][8] != data[1][8]:
			if begin_r - end_r < thr:
				classi = "DEL_INV"
				dict_number[classi] += 1
			elif begin_r - end_r >= thr:
				classi = "INS_INV"
				dict_number[classi] += 1
			
		for cosa in data:
			c += 1
			match_reads[el][c] = [classi] + cosa + ["DEF"]
			cosa = "\t".join(match_reads[el][c])
			outnf.write(cosa + "\n")

out.close()
outnf.close()

#print match_reads

# calculating weighted average of divergence
dict_c = {}

for lis in match_reads:
	for li in match_reads[lis]:
		clas = li[0] #class
		div = float(li[2])
		begin = int(li[6])
		end = int(li[7])+1
		f = li[-1]
		li_info = [div, end-begin, f]

		# filling dict_c with divergence and length data
		look = clas in dict_c
		if look == True:
			dict_c[clas].append(li_info)

		elif look == False:
			dict_c[clas] = [li_info]

func = ["FULL", "EXT5", "EXT3"]

numerator_all = 0
numerator_func = 0
numerator_def = 0
divisor_all = 0
divisor_func = 0
divisor_def = 0

#print dict_c

for eli in dict_c:
	el = dict_c[eli] # for each class in eli
	numerator = 0
	divisor = 0
	for e in el:
		print e
		numerator += e[0]*e[1]
		numerator_all += e[0]*e[1]
		divisor += e[1]
		divisor_all += e[1]
#		if eli in func: #If FULL, EXT3 or EXT5
		if e[-1] == "FUN": #If FULL, EXT3 or EXT5
			numerator_func += e[0]*e[1]
			divisor_func += e[1]
		else:
			numerator_def += e[0]*e[1]
			divisor_def += e[1]

	pond = 1.0 * numerator/divisor
	pond = str(round(pond,1))
	dict_pond[eli] = pond

try:
	a = dict_c[""]
	print "Sequeneces do not classified"
except:
	pass

pond_all = 1.0 * numerator_all/divisor_all
pond_all = str(round(pond_all,1))
dict_pond["ALL"] = pond_all

pond_func = 1.0 * numerator_func/divisor_func
pond_def = 1.0 * numerator_def/divisor_def

print "FUNC_NUC: %s\nDEF_NUC: %s" % (str(divisor_func), str(divisor_def))
print "ALL_NUC: %s" % str(nuc)
print "FUNC: %s\nDEF: %s" % (str(round(pond_func,1)), str(round(pond_def,1))) 

# Get sums for all sequences
number_all = 0
trans_all = 0
trans_all_nuc = 0
for el in dict_number:
	number_all += dict_number[el]
	try:
		trans_all += dict_trans[el]
		trans_all_nuc += dict_trans_nuc[el]
	except:
		pass
dict_number["ALL"] = number_all
dict_trans["ALL"] = trans_all
dict_trans_nuc["ALL"] = trans_all_nuc

print dict_pond
print dict_number
print dict_trans
print dict_trans_nuc

# creating summary file
sum_file = open(data_file + ".sum", "w")

sum_file.write("Name\tALL\tFULL\tEXT5\tEXT3\tDEL\tINS\tDEL_INV\tINS_INV\n")
sum_file.write(data_file)

classes = ["ALL", "FULL", "EXT5", "EXT3", "DEL", "INS", "DEL_INV", "INS_INV"]

for cl in classes:
	result = ""
	result = result + str(dict_number[cl])
	result = result + " (" + str(dict_pond[cl]) + ")"
	try:
		result = result + " " + str(dict_trans[cl])
	except:
		pass
	sum_file.write("\t" + result)

sum_file.write("\n")

sum_file.close()


print "TRUNC_5P: " + str(number_trunc_5p) + " " + str(number_trunc_5pp) + " " + str(number_trunc_5p_fun)
print "TRUNC_3P: " + str(number_trunc_3p) + " " + str(number_trunc_3pp) + " " + str(number_trunc_3p_fun)

#####################################

for line in file[3:]:
	text = line.split()
	begin = int(text[5])-1
	end = int(text[6])
	left = text[7]
	left = left.split("(",1)[1].split(")")[0]
	left = int(left)
#	if end - begin > 96:
	if begin + left < 6:
		diver = float(text[1])
		length = end-begin
		annot = text[9]
		look = annot in dict_count
		if look == True:
			dict_count[annot] += 1
			dict_diver[annot].append(diver)
			dict_allen[annot].append(length)
		elif look == False:
			dict_count[annot] = 1
			dict_diver[annot] = [diver]
			dict_allen[annot] = [length]

print "ELI FUN: %s\nELI DEF: %s\n" % (str(eli_fun), str(eli_def))

print dict_count


