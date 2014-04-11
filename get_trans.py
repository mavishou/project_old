#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

'''
This script is used to get the entire transcript annotation from the V4expCaculate_combined_houm.gtf file
'''

import modify_gtf as gtf
import sys

def summaryLast():
	global curTransId
	global curGeneId
	global curCh
	global curSrc
	global curStrand
	global lCoords

	# pool all coordinates together, the min is the first, the max is the last
	firstCoord = min(lCoords)
	lastCoord = max(lCoords)

	lOut = [curCh, curSrc, 'transcript', firstCoord, lastCoord, '.', '-', '.', 'gene_id ' + curGeneId + '; transcript_id ' + curTransId]
	print '\t'.join(lOut)


curTransId = ''

for line in sys.stdin.readlines():
	li = line.rstrip('\n').split('\t')
	dAnnos = gtf.processAnnotation(li[8])
	transId = dAnnos['transcript_id']
	geneId = dAnnos['gene_id']

	if transId != curTransId and curTransId != '':
		summaryLast()

	if transId != curTransId:
		curTransId = transId
		curGeneId = geneId
		curCh, curSrc = li[0:2]
		curStrand = li[6]
		lCoords = []

	lCoords += li[3:5]

summaryLast()



