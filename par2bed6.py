# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
# remove the first line
firstLine = sys.stdin.readline()

for line in sys.stdin.readlines():
	line = line.rstrip('\n')
	li = line.split('\t')
	start = str(int(li[2])-1)
	lOutput = [li[0], start, li[3], li[4], '.', li[1]]
	sOutput = '\t'.join(lOutput)
	print sOutput


