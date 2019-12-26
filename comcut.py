#!/usr/bin/python

import sys,re,pickle,getopt,urllib2,os,subprocess
import numpy as np
import time

usage = """

*****************************************************************************

This script.....

Usage

	python cCUT.py <-option argument> <amino acid sequence or Fasta file>


Options:

	-s (--species) <Species>	Comma-sep species names (default: Homo_sapiens)
	-w (--weights) <n>		Comma-sep integers/floats (default: 1)
	-m (--maxhits) <n>		Max hits per input sequence (default: 5)

	-r (--removerare)		Remove codons with frequency < 0.05 (default: False)
	-c (--reduceCG)			Reduce occurrence of CpG repeats (default: False)
	-e (--avoidEnzymes) <Enzymes>	Comma-sep enzyme names to avoid (default: None)
	-d (--reduceDonor)		Reduce occurrence of 5' donor splice sites (default: False)
	-a (--reduceAcceptor)		Reduce occurrence of 3' acceptor splice sites (default: False)	

	-f (--fullanalysis)		Save list of hits (default: False)
	-o (--output)			Save compromised codon table (default: False)
	-h (--help)			Display script usage

*****************************************************************************

"""

########################################
############ Class and Def #############
########################################

class Params():
	def __init__(self):

		# define necessary variables
		self.list_organisms = ['Homo_sapiens']
		self.list_weights = [1]
		self.bool_savetable = False
		self.bool_removerare = False
		self.bool_reduceCG = False
		self.list_RE = []
		self.bool_check5prime = False
		self.bool_check3prime = False
		self.hits = 5
		self.fullanalyze = False

		# parse command-line arguments and check for proper input
		self.argparse()

	def argparse(self):

		# attempt to read command-line option-argument tuples and mandatory argument.
		try:
			options, rawinput = getopt.getopt(sys.argv[1:], "m:s:w:orce:adfh",["maxhits","species=","weights=","output","removerare","reduceCG","avoidEnzymes","reduceDonor","reduceAcceptor","fullanalysis","help"])
			self.rawinput = str(rawinput[0]).strip()
		except getopt.GetoptError, err:
			sys.exit(usage)		# exit script and print usage if arguments to options are not provided 
		except IndexError:
			sys.exit(usage)		# exit script and print usage if command-line input is not provided

		# parse command-line options into its appropriate variables/actions
		for opt,arg in options:
			if opt in ("-s","--species"):
				self.list_organisms = map(str,arg.split(','))
				self.list_weights = [1]*len(self.list_organisms)
			if opt in ("-w","--weights"):
				self.list_weights = map(float,arg.split(','))
			if opt in ("-o","--output"):
				self.bool_savetable = True
			if opt in ("-r","--removerare"):
				self.bool_removerare = True
			if opt in ("-c","--reduceCG"):
				self.bool_reduceCG = True
			if opt in ("-e","--avoidEnzymes"):
				self.list_RE = map(str,arg.split(','))
			if opt in ("-a","--reduceDonor"):
				self.bool_check5prime = True
			if opt in ("-d","--reduceAcceptor"):
				self.bool_check3prime = True
			if opt in ("-m","--maxhits"):
				self.hits = int(arg)
			if opt in ("-f","--fullanalysis"):
				self.fullanalyze = True
			if opt in ("-h","--help"):
				sys.exit(usage)

def argcheck():
	# Additional preparations and checks
	dict_codonrank,aa = pickle.load(open('matrix','r'))		#pre-load dictionaries and aa list
	dict_species = pickle.load(open('species','r'))
	dict_RE = pickle.load(open('RE_matrix','r'))
	if not re.search(r'.txt',Parameters.rawinput):
		if  re.search(r'[^GALMFWKQESPVICYHRNDT*]',Parameters.rawinput.upper()): sys.exit('Error: Input sequence may contain illegal characters')			# check for non-amino acid characters
	if len(Parameters.list_weights) != len(Parameters.list_organisms): sys.exit('Error: weights do not match number of organisms')			# check for correct number of weights
	
	return dict_codonrank,dict_species,dict_RE

def prepanalysis():
	list_spnumber = preporganism(Parameters.list_organisms)		# provide organismID for species without locally-saved table
	prepmatrix(list_spnumber)
	processmatrix()

def preporganism(organisms):
	list_speciesID = []
	prelistorganisms = os.listdir('/Users/Fursham/Documents/Programming/Project/010_CCUT/Codon_tables')
	for eachorganism in organisms:
		if eachorganism in prelistorganisms:
			list_speciesID.append(eachorganism)
		else:
			list_speciesID.append(dict_species[re.sub('_', ' ',eachorganism)])

	return list_speciesID

def prepmatrix(orgID):
	for thisorgID in orgID:			# iterate list of organisms/organismID

		# retrieve codon table for selected organisms, either from locally-saved tables or from kazusa.or.jp
		if re.match(r'[\d]+',thisorgID):		# from kazusa
			url = "".join(["http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=",thisorgID,"&aa=1&style=N"])
			response = urllib2.urlopen(url)
			source = response.read()
			list_codonusage = list(re.split('\n|\)  ',(re.split('\n<PRE>\n|\n</PRE>\n',source)[1])))
		else:			# from local tables
			fileoftable = open('Codon_tables/%s'%thisorgID,'r')
			table = fileoftable.read()
			list_codonusage = list(re.split('\n|\)  ',table))

		# this loop will extract the frequency usage of each codons and save it to dict_codonrank
		for eachline in list_codonusage:
			if eachline == '':continue
			list_data = list(re.split(' |  ',eachline))
			list_data = filter(None,list_data)
			dict_codonrank[list_data[1]][list_data[0]].append([float(list_data[3])* Parameters.list_weights[orgID.index(thisorgID)],(float(list_data[3])/float(list_data[2]))])

def processmatrix():

	# this loop will rank the codons based on frequency of usage in this organism, apply weights to it and sum it to the dict_codonrank
	# each codon is a key to a list with structure [[occurrence_per_thousand,totalnumberforAA]*eachorganisms,mean_of_weighted_frequency,probability,probability_cumulative]
	for eachAA in dict_codonrank:
		sumofoccurrence = 0		
		baseprobability = 0

		# loop to calculate mean of occurrence, append it to the list and add it to sumofoccurrence for each AA
		for eachCodon in dict_codonrank[eachAA]:
			thismean = np.mean(map(lambda x:dict_codonrank[eachAA][eachCodon][x][0],xrange(len(dict_codonrank[eachAA][eachCodon]))))
			dict_codonrank[eachAA][eachCodon].append(thismean)
			sumofoccurrence += thismean

		# loop to calculate new probability and append it to the list
		for eachCodon in dict_codonrank[eachAA]:
			dict_codonrank[eachAA][eachCodon].append((dict_codonrank[eachAA][eachCodon][-1]/sumofoccurrence))

			# remove rare codons if requested, and calculate cumulative probability
			if Parameters.bool_removerare == False or any(i > 0.1 for i in (filter(None,map(lambda x: (x[0]/x[1]) if type(x) is list else '',dict_codonrank[eachAA][eachCodon])))):
				dict_codonrank[eachAA][eachCodon].append((dict_codonrank[eachAA][eachCodon][-2]/sumofoccurrence)+baseprobability)
				baseprobability += dict_codonrank[eachAA][eachCodon][-3]/sumofoccurrence
			elif Parameters.bool_removerare == True and any(i < 0.1 for i in (filter(None,map(lambda x: (x[0]/x[1]) if type(x) is list else '',dict_codonrank[eachAA][eachCodon])))):

				dict_codonrank[eachAA][eachCodon].append(0)
			
	if Parameters.bool_savetable == True: savematrix();		# output the compromised matrix if requested

def savematrix():		
	output_handle = open('_'.join(['compromisedcodon','_'.join(Parameters.list_organisms)]),'w')		# prepare file for writing

	# prints out details of the new compromised table and important headers
	output_handle.write("Organisms used to generate compromised codon usage: %s\n" %Parameters.list_organisms)
	output_handle.write("Weightage: %s\n\n" %Parameters.list_weights)
	output_handle.write("Amino Acid\t%s\tCompromised frequency\n" %('\t'.join(Parameters.list_organisms)))

	# loop through dictionary and prints out the sum of weighted ranks for each codon
	for thisAA in dict_codonrank:
		for thisCodon in dict_codonrank[thisAA]:
			initialfreq = '\t'.join(filter(None,map(lambda x: str(x[0]/x[1]) if type(x) is list else '', dict_codonrank[thisAA][thisCodon])))

			output_handle.write("%s\t%s\t" %(thisAA,thisCodon))
			output_handle.write("%s\t" %(initialfreq))
			output_handle.write("%.2f\t"% (dict_codonrank[thisAA][thisCodon][-2]))	
			output_handle.write("%.3f\n" %dict_codonrank[thisAA][thisCodon][-1])
		output_handle.write("\n")

def processinput():
	# to process fasta file and output optimized codon into a new fasta file
	if re.search(r'.txt',Parameters.rawinput):
		outputfile = "".join([Parameters.rawinput[:-4],"_codonoptimized.txt"])		# sweet way to append existing .txt files

		# prepare files for reading and writing
		try:input_handle = open(Parameters.rawinput,"r")
		except IOError:sys.exit("Error. Input file might not be there!")
		output_handle = open (outputfile,"w")
		inputseq = ''

		# loop to read fasta input lines
		for eachline in input_handle:
			eachline = eachline.rstrip()		# remove newlines

			if re.match(r'^>',eachline):	# check for fasta handlename
				if len(inputseq) > 0:
					candidateseq = {}
					codonseq = choosecodons(inputseq)
					response,CGratio,REfound,list5ss,list3ss = checksequence(codonseq)
					if response == 'pass':
						candidateseq[codonseq] = {}
						candidateseq[codonseq]['CG'] = CGratio
						candidateseq[codonseq]['RE'] = REfound
						candidateseq[codonseq]['5ss'] = list5ss[:]
						candidateseq[codonseq]['3ss'] = list3ss[:]
					bestseq = ''.join([y if ((candidateseq[y]['5ss'][0] + candidateseq[y]['3ss'][0]) == min([(candidateseq[x]['5ss'][0] + candidateseq[x]['3ss'][0]) for x in candidateseq]) )else '' for y in candidateseq])
					output_handle.write("\nBest codon sequence:\n%s\nCG ratio: %.3f\tRE found: %s\tMax 5'ss: %s\tMax 3'ss: %s\n" %(bestseq,candidateseq[bestseq]['CG'],candidateseq[bestseq]['RE'],max(candidateseq[bestseq]['5ss'][1]),max(candidateseq[bestseq]['3ss'][1])))

					inputseq = ''
					output_handle.write("%s\n" %eachline.replace(">",">codopt_"))	# pop the string "codopt_" in the fasta handlename and print onto output_handle


				else:
					output_handle.write("%s\n" %eachline.replace(">",">codopt_"))	# pop the string "codopt_" in the fasta handlename and print onto output_handle

			else:		# check if line contain illegal characters and if no, output the optimized codon sequence
				if  re.search(r'[^GALMFWKQESPVICYHRNDT*]',eachline.upper()): sys.exit('Error: Fasta file may contain illegal characters')		
				inputseq += eachline.upper()

	# to process sequences from command-line and output into command
	else: 
		candidateseq = {}
		while len(candidateseq) <= Parameters.hits:
		# for i in xrange(1):
			thisseq = choosecodons(Parameters.rawinput.upper())
			response,CGratio,REfound,list5ss,list3ss = checksequence(thisseq)
			# print response
			if response == 'pass':
				candidateseq[thisseq] = {}
				candidateseq[thisseq]['CG'] = CGratio
				candidateseq[thisseq]['RE'] = REfound
				candidateseq[thisseq]['5ss'] = list5ss[:]
				candidateseq[thisseq]['3ss'] = list3ss[:]
		prepoutput(candidateseq)

def choosecodons(seq):
	return re.sub(r'U','T', ''.join(map(lambda x: next(y for y in dict_codonrank[x] if dict_codonrank[x][y][-1]>(np.random.uniform(0,max(map(lambda z:dict_codonrank[x][z][-1],dict_codonrank[x]))))),list(seq))))		#prettyyyyyy

def checksequence(newseq):

	CGratio = checkCG(newseq)
	REfound = checkRE(newseq)
	list5ss = check5(newseq)
	list3ss = check3(newseq)

	# print CGratio,REfound,(list5ss),(list3ss)

	if (CGratio < 10) and (REfound == 0):
		if (len(list5ss[1]) == 0) and (len(list3ss[1]) == 0):
			return 'pass',CGratio,REfound,list5ss,list3ss
		elif (max(list5ss[1]) < 8 and max(list3ss[1]) < 8):
			return 'pass',CGratio,REfound,list5ss,list3ss
		else: return 'fail',0,0,0,0
	else: return 'fail',0,0,0,0

def checkCG(newseq):

	if (Parameters.bool_reduceCG):
		CGrepeat = (len(re.findall('CG',newseq))/float(len(newseq)))*100
		# CGrepeat = len(re.findall('CG',newseq))
		CGtotal = len(re.findall('C',newseq))/float(len(newseq)) * len(re.findall('G',newseq))/float(len(newseq))
		return CGrepeat
	else: return 0

def checkRE(newseq):
	RE = 0
	for eachRE in Parameters.list_RE:
		thisRE = [key for key in dict_RE if re.search(eachRE,key)]
		REseq = prepRE(dict_RE[''.join(thisRE)])
		RE += len(re.findall(REseq,newseq))
	return RE

def check5(newseq):

	if (Parameters.bool_check5prime):
		list_seq = map(lambda x: ''.join(newseq[x:x+9]),xrange(0,len(newseq)-9))
		donorcout = subprocess.Popen(['perl','score5.pl',','.join(list_seq)],stdout=subprocess.PIPE)
		list_donorvalues,value = donorcout.communicate()
		list_donorvalues = map(float,filter(None,list_donorvalues.split('\t')))

		list_fivess = filter(None,(map(lambda x: x if x > 3 else 0.0001,list_donorvalues)))
		list_fivessseq = filter(None,(map(lambda x: ''.join(list(newseq)[x:x+9]) if list_donorvalues[x] > 3 else '',xrange(0,len(list_donorvalues)-9))))
		list_fivesspos = filter(None,(map(lambda x: x+1 if list_donorvalues[x] > 3 else '',xrange(0,len(list_donorvalues)-9))))

		return [sum(list_fivess),list_fivess,list_fivessseq,list_fivesspos]
	else: return [0,[0],[0],[0]]

def check3(newseq):

	if (Parameters.bool_check3prime):
		listofseq = map(lambda x: ''.join(newseq[x:x+23]),xrange(0,len(newseq)-23))
		acceptor = subprocess.Popen(['perl','score3.pl',','.join(listofseq)],stdout=subprocess.PIPE)
		list_acceptorvalues,value = acceptor.communicate()
		list_acceptorvalues = map(float,filter(None,list_acceptorvalues.split('\t')))

		list_threess = filter(None,map(lambda x: x if x > 3 else 0.0001,list_acceptorvalues))
		list_threessseq = filter(None,(map(lambda x: ''.join(list(newseq)[x:x+23]) if list_acceptorvalues[x] > 3 else '',xrange(0,len(list_acceptorvalues)-23))))
		list_threesspos = filter(None,(map(lambda x: x+1 if list_acceptorvalues[x] > 3 else '',xrange(0,len(list_acceptorvalues)-23))))

		return [sum(list_threess),list_threess,list_threessseq,list_threesspos]
	else: return [0,[0],[0],[0]]

def prepRE(seq):
	seq = re.sub('B','[CGT]',seq)
	seq = re.sub('D','[AGT]',seq)
	seq = re.sub('H','[ACT]',seq)
	seq = re.sub('K','[GT]',seq)
	seq = re.sub('M','[AC]',seq)
	seq = re.sub('N','[ACGT]',seq)
	seq = re.sub('R','[AG]',seq)
	seq = re.sub('S','[CG]',seq)
	seq = re.sub('V','[ACG]',seq)
	seq = re.sub('W','[AT]',seq)
	seq = re.sub('Y','[CT]',seq)
	return seq

def prepoutput(candidateseq):
	bestseq = np.random.choice(filter(None,([y if ((candidateseq[y]['5ss'][0] + candidateseq[y]['3ss'][0] + candidateseq[y]['CG']) == min([(candidateseq[x]['5ss'][0] + candidateseq[x]['3ss'][0] + candidateseq[x]['CG']) for x in candidateseq])) else '' for y in candidateseq])))
	
	print "\nBest codon sequence:\n%s\nCG ratio: %.3f\tRE found: %s\tMax 5'ss: %s\tMax 3'ss: %s\n" %(bestseq,candidateseq[bestseq]['CG'],candidateseq[bestseq]['RE'],max(candidateseq[bestseq]['5ss'][1]),max(candidateseq[bestseq]['3ss'][1]))

	if Parameters.fullanalyze == True:

		out_handle = open('cCUT_output.txt','w')
		out_handle.write("Organisms used to generate compromised codon usage: %s\n" %Parameters.list_organisms)
		out_handle.write("Weightage: %s\n\n" %Parameters.list_weights)

		out_handle.write("Best codon sequence:\n%s\nCG ratio: %.3f\tRE found: %s\tMax 5'ss: %s\tMax 3'ss: %s\n\n" %(bestseq,candidateseq[bestseq]['CG'],candidateseq[bestseq]['RE'],(candidateseq[bestseq]['5ss'][0]),(candidateseq[bestseq]['3ss'][0])))
		out_handle.write("List of all generated codon sequence\n")
		for eachkey in candidateseq:
			out_handle.write("Sequence:\n%s\nCG ratio: %.3f\tRE found: %s\tMax 5'ss: %s\tMax 3'ss: %s\n\n" %(eachkey,candidateseq[eachkey]['CG'],candidateseq[eachkey]['RE'],(candidateseq[eachkey]['5ss'][0]),(candidateseq[eachkey]['3ss'][0])))


########################################
################ Main ##################
########################################

if __name__ == "__main__":
	start_time=time.time()

	Parameters = Params()
	dict_codonrank,dict_species,dict_RE = argcheck()	# parse command-line arguments and define global variables
	prepanalysis(); candidateseq = processinput()		# process the input seq/file and output codon-optimized DNA sequence

	elapsed_time=time.time()-start_time
	print "Time taken for analysis: %.3f seconds" %elapsed_time



