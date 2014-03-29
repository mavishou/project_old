# -*- coding: utf-8 -*
#!/usr/bin/python

import sys
import numpy as np
import pdb
import cPickle
import copy

class Transcript:
	#pdb.set_trace()
	def __init__(self, geneId, transID, curStrand, exons, cdss):
		self.geneId = geneId
		self.transID = transID
		self.exons = exons
		self.cdss = cdss
		self.strand = curStrand
		self.generalInfo = self.exons[0][0:2] + self.exons[0][5:8]
		self.exons, self.exonCods = self.getCoordiates(self.exons)
		self.outputList = []

		self.outputList.extend(self.__getExonOut())
		# self.exons[0].append('test')
# 		# exon大于1才会有intron
# 		if len(self.exons) > 1:
# 			self.intronCods = self.getIntronCords(self)
# 			self.outputList.extend(self.getIntronOut(self))

# 		# 如果有cds
# 		if len(self.cdss) > 0:
# 			self.cdsCods = self.getCoordiates(self, self.cdss)
# 			self.utr5Cods, self.utr3Cods = self.getUTR(self)
# 			self.outputList.extend(getCDSOut(self))
# 			self.outputList.extend(getUTROut(self))
# 		self.getFinalOut(self)

	def getCoordiates(self, originList):
		coordinates = []
		for l in originList:
			coordinates.append([int(s) for s in l[3:5]])
		coordinates = np.array(coordinates)
		idx = np.lexsort((coordinates[:,1],coordinates[:,0]))
		coordinates = coordinates[idx]
		# orinigal list同时进行排序，否则和坐标不同步了
		# 注意改变这里并不会改变
		# originList = [originList[i] for i in idx]
		originList = copy.deepcopy([originList[i] for i in idx])
		return originList, coordinates

# 	def getIntronCords(self):
# 		'''input: exonCods'''
# 		intronCods = self.exonCods.reshape(1,-1)[:, 1:-1].reshape(-1,2)
# 		intronCods[:, 0] = intronCods[:, 0] + 1
# 		intronCods[:, 1] = intronCods[:, 1] - 1
# 		return intronCods

	
# 	def cmpExonCDS(self, exon, cds, mark):
# 		'''input: 
# 				a pair of exon coordinates
# 				a pair of cds coordinates
# 				mark: 0 -> begin; 1 -> end
# 		'''
# 		if (mark == 0 and exon[0] < cds[0] and exon[1] == cds[1]) or (mark == 1 and exon[1] > cds[1] and exon[0] == cds[0]):
# 			return True
# 		else:
# 			return False

# 	def getUTR(self):
# 		'''input:exonCods, cdsCods
# 		'''
# 		utr3 = []
# 		utr5 = []
# 		cdscoords = cdsCods.copy().reshape(1,-1)[0,]
# 		cdscoords.sort()
# 		cdsFirst = cdscoords[0] - 1
# 		cdsLast = cdscoords[-1] + 1
# 		exoncoords = exonCods.copy().reshape(1,-1)[0,]
# 		exoncoords.sort()
# 		exonFirst = exoncoords[0]
# 		exonLast = exoncoords[-1]
# 		if cdsFirst < exonFirst or cdsLast > exonLast:
# 			sys.stderr.write('The CDS is beyond the exon! ' + self.transID + '\n')
# 			sys.exit(1)
# 		for ec in exonCods:
# 			bf1 = getOverlap(ec, [exonFirst, cdsFirst]):
# 			if bf1 != False:
# 				utr5.append(bf1)
# 			bf2 = getOverlap(ec, [cdsLast, exonLast])
# 			if bf2 != False:
# 				utr3.append(bf2)
# 		utr3 = np.array(utr3)
# 		utr5 = np.array(utr5)
# 		if self.strand == '-':
# 			utr5, utr3 = utr3, utr5
# 		return utr5, utr3


# 	def getUTROut(self):
# 		'''input: utrCods, strand'''
# 		annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID
# 		utrOut=[]
# 		for utr in self.utr5Cods:
# 			utrOut.append(self.generalInfo[:2] + "5_utr" + [str(s) for s in utr] + self.generalInfo[2:] + annotation)
# 		for utr in self.utr3Cods:
# 			utrOut.append(self.generalInfo[:2] + "3_utr" + [str(s) for s in utr] + self.generalInfo[2:] + annotation)
# 		return(utrOut)	

	def __getExonOut(self):
		'''input: exons, stand'''
		# 如果是负链，exon颠倒
		if self.strand == '-':
			# exonsOut = self.exons[::-1]
			exonsOut = copy.deepcopy(self.exons)[::-1]
		else:
			exonsOut = copy.deepcopy(self.exons)
		for i in range(len(exonsOut)):
			exonNum = i + 1
			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; exon_number ' + str(exonNum)
			exonsOut[i].append(annotation)
		return(exonsOut)

# 	def getCDSOut(self):
# 		'''input: cdss, stand'''
# 		if self.strand == '-':
# 			cdssOut = self.cdss[::-1]
# 		else:
# 			cdssOut = self.cdss[:]
# 		for i in range(len(cdssOut)):
# 			cdsNum = i + 1
# 			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; CDS_number ' + cdsNum
# 			cdssOut[i].append(annotation)
# 		return(getExonOut)

# 	def getIntronOut(self):
# 		'''input: intronCods, info, stand'''
# 		if self.strand == '-':
# 			myIntronCods = self.intronCods[::-1]
# 		else:
# 			myIntronCods = self.intronCods[:]
# 		intronOut = []
# 		for i in range(len(myIntronCods)):
# 			intronNum = i + 1
# 			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; intron_number ' + intronNum
# 			out = self.generalInfo[:2] + "intron" + [str(s) for s in myIntronCods[i]] + self.generalInfo[2:] + annotation
# 			intronOut.append(out)
# 		return(intronOut)

	
		
# 	def getFinalOut(self):
# 		'''input: outputList, strand'''
# 		tmp = self.getCoordiates(self, self.outputList)
# 		if self.strand == '-':
# 			self.outputList = self.outputList[::-1]

curtTransId = ''
curGeneId = ''
curStrand = ''
exons = []
cdss = []

def getOverlap(x, y):
	if x[1] > y[0] or y[1] > x[0]:
		return False
	else:
		coordinates = np.append(x, y)
		return coordinates.sort()[1:3]

def processAnnotation(annotation):
	'''input: the 9th colume of a gtf file
	output: a dict of annotations 
	'''
	ants = annotation.split(';')
	# 如果以分号结尾，则最后一个元素为空，把空的去掉，并把每个元素开头的空格去掉
	ants = [a.strip() for a in ants if a != '']
	antsDict = {}
	for a in ants:
		antSplit=a.split(' ')
		##check length
		antsDict[antSplit[0]] = antSplit[1]
	return antsDict

def finalOutput(myList):
	for l in myList:
		outputLine = '\t'.join(l)
		print outputLine

####################################################
for line in sys.stdin.readlines():
	# pdb.set_trace()
	li = line.rstrip('\n').split('\t')
	# annotation
	antsDict = processAnnotation(li[8])
	
	# check whether have gene id and trans id
	if not 'gene_id' in antsDict:
		sys.stderr.write('No gene_id!\n')
		print line
		sys.exit(1)
	if not 'transcript_id' in antsDict:
		sys.stderr.write('No transcript_id!\n')
		print line
		sys.exit(1)

	if antsDict['transcript_id'] != curtTransId and curtTransId !='': ##### 1st
		trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
		finalOutput(trans.outputList)

	if antsDict['transcript_id'] != curtTransId: ##### 2nd
		curtTransId = antsDict['transcript_id']
		curGeneId = antsDict['gene_id']
		curStrand = li[6]
		exons=[]
		cdss=[]
	
	if li[2] == 'exon':
		exons.append(li[:8])
	if li[2] == 'CDS':
		cdss.append(li[:8])

##test beign
# print str(exons)
# print str(cdss)
# print curStrand
# print curtTransId
# print curGeneId
##test end

##test beign
# w = open('exons.tmp', 'w')
trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
# cPickle.dump(trans.exons, w)
# print len(li)
# print li[:8]
print exons
print
print trans.exons
print
print trans.exonCods
print
print trans.outputList
print
# print trans.cdss
print trans.transID
print trans.geneId
# print trans.generalInfo
##test end


# the last output
# trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
# finalOutput(trans.outputList)



