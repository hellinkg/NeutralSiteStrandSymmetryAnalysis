'''
CountRepeatType.py
@author: Kristin Helling
'''
from RMfileReader import *
import numpy
		
def createMatrix(queryMem, TorFtrans, TorFCG):
	""" 
	Method that takes in a Repeat object, aligns the two sequences, and
	returns the dimension of the matrix and a transition matrix as a string.
	queryMem: Repeat object created by RMfileReader.py.
	TorFtrans: Boolean value for only counting transitions in the matrix.
	TorFCG: Boolean value for only counting non-CpG sites.
	"""	

	CorPseq = queryMem.CorPseq  # ancestral sequence
	seq = queryMem.seq  # modern sequence
	
	# Matrix H = ACGT == 0,1,2,3 modern sequence
	# Matrix V = ACGT == 0,1,2,3 ancestor sequence

	countMat = numpy.zeros((4, 4), dtype=numpy.int)

	count = 0
	limit = 0
	if len(seq) <= len(CorPseq):
		limit = len(seq)
	else:
		limit = len(CorPseq)

	while(count < limit):
		#!!!!!!!!!!!! only does C or + difference for False and False instance currently!!!!!!!!!!
		if TorFtrans == True and TorFCG == True:
			if CorPseq[count] == "A" and seq[count] == "G":
				countMat[0][2] += 1
			elif CorPseq[count] == "C":
				if (count > 0 and count < limit - 1) and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G") and seq[count] == "T":
					countMat[1][3] += 1  
				elif count == 0 and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G") and seq[count] == "T":
					countMat[1][3] += 1	
				elif count == limit - 1 and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and seq[count] == "T":
					countMat[1][3] += 1			
				else:
					if CorPseq[count] == "T":
						countMat[1][3] += 1   	   
			elif CorPseq[count] == "G":
				if (count > 0 and count < limit - 1) and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G") and seq[count] == "A":
					countMat[2][0] += 1
				elif count == 0 and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G") and seq[count] == "A":
					countMat[2][0] += 1
				elif count == limit - 1 and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and seq[count] == "A":
					countMat[2][0] += 1
				else:
					if seq[count] == "A":
						countMat[2][0] += 1
			elif CorPseq[count] == "T" and seq[count] == "C":
				countMat[3][1] += 1

		elif TorFtrans == False and TorFCG == True:
			if CorPseq[count] == "A":
				if seq[count] == "A":
					countMat[0][0] += 1
				elif seq[count] == "C":
					countMat[0][1] += 1
				elif seq[count] == "G":
					countMat[0][2] += 1
				elif seq[count] == "T":
					countMat[0][3] += 1
			
			elif CorPseq[count] == "C":
				if (count > 0 and count < limit - 1) and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G"):
					if seq[count] == "A":
						countMat[1][0] += 1
					elif seq[count] == "C":
						countMat[1][1] += 1
					elif seq[count] == "G":
						countMat[1][2] += 1
					elif seq[count] == "T":
						countMat[1][3] += 1  
				elif count == 0 and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G"):
					if seq[count] == "A":
						countMat[1][0] += 1
					elif seq[count] == "C":
						countMat[1][1] += 1
					elif seq[count] == "G":
						countMat[1][2] += 1
					elif seq[count] == "T":
						countMat[1][3] += 1	
				elif count == limit - 1 and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G"):
					if seq[count] == "A":
						countMat[1][0] += 1
					elif seq[count] == "C":
						countMat[1][1] += 1
					elif seq[count] == "G":
						countMat[1][2] += 1
					elif seq[count] == "T":
						countMat[1][3] += 1			
				else:
					if seq[count] == "A":
						countMat[1][0] += 1
					elif seq[count] == "C":
						countMat[1][1] += 1
					elif seq[count] == "G":
						countMat[1][2] += 1
					elif seq[count] == "T":
						countMat[1][3] += 1   
				
			elif CorPseq[count] == "G":
				if (count > 0 and count < limit - 1) and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G") and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G"):
					if seq[count] == "A":
						countMat[2][0] += 1
					elif seq[count] == "C":
						countMat[2][1] += 1
					elif seq[count] == "G":
						countMat[2][2] += 1
					elif seq[count] == "T":
						countMat[2][3] += 1
				elif count == 0 and (CorPseq[count + 1] != "C" or CorPseq[count + 1] != "G"):
					if seq[count] == "A":
						countMat[2][0] += 1
					elif seq[count] == "C":
						countMat[2][1] += 1
					elif seq[count] == "G":
						countMat[2][2] += 1
					elif seq[count] == "T":
						countMat[2][3] += 1	
				elif count == limit - 1 and (CorPseq[count - 1] != "C" or CorPseq[count - 1] != "G"):
					if seq[count] == "A":
						countMat[2][0] += 1
					elif seq[count] == "C":
						countMat[2][1] += 1
					elif seq[count] == "G":
						countMat[2][2] += 1
					elif seq[count] == "T":
						countMat[2][3] += 1	
				elif CorPseq[count] == "T":
					if seq[count] == "A":
						countMat[3][0] += 1
					elif seq[count] == "C":
						countMat[3][1] += 1
					elif seq[count] == "G":
						countMat[3][2] += 1
					elif seq[count] == "T":
						countMat[3][3] += 1	
				else:
					if seq[count] == "A":
						countMat[2][0] += 1
					elif seq[count] == "C":
						countMat[2][1] += 1
					elif seq[count] == "G":
						countMat[2][2] += 1
					elif seq[count] == "T":
						countMat[2][3] += 1		

		elif TorFtrans == True and TorFCG == False:
			if CorPseq[count] == "A" and seq[count] == "G":
				countMat[0][2] += 1
			elif CorPseq[count] == "C" and seq[count] == "T":
				countMat[1][3] += 1   	   
			elif CorPseq[count] == "G" and seq[count] == "A":
				countMat[2][0] += 1
			elif CorPseq[count] == "T" and seq[count] == "C":
				countMat[3][1] += 1	
		
		else:
			if CorPseq[count] == "A":
				if seq[count] == "A":
					countMat[0][0] += 1
				elif seq[count] == "C":
					countMat[0][1] += 1
				elif seq[count] == "G":
					countMat[0][2] += 1
				elif seq[count] == "T":
					countMat[0][3] += 1
				
			elif CorPseq[count] == "C":
				if seq[count] == "A":
					countMat[1][0] += 1
				elif seq[count] == "C":
					countMat[1][1] += 1
				elif seq[count] == "G":
					countMat[1][2] += 1
				elif seq[count] == "T":
					countMat[1][3] += 1   
				
			elif CorPseq[count] == "G":
				if seq[count] == "A":
					countMat[2][0] += 1
				elif seq[count] == "C":
					countMat[2][1] += 1
				elif seq[count] == "G":
					countMat[2][2] += 1
				elif seq[count] == "T":
					countMat[2][3] += 1		
				
			elif CorPseq[count] == "T":
				if seq[count] == "A":
					countMat[3][0] += 1
				elif seq[count] == "C":
					countMat[3][1] += 1
				elif seq[count] == "G":
					countMat[3][2] += 1
				elif seq[count] == "T":
					countMat[3][3] += 1	  	
	
		count += 1
		
	outLine = str(len(countMat)) + " "   
		
	for y in range(4):
		for z in range(4):
			outLine = outLine + str(countMat[y][z]) + " "
	
	return outLine