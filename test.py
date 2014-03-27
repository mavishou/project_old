#!/usr/bin/python
import sys
for line in sys.stdin.readlines():
	line = line.rstrip('\n')
	li = line.split('\t')
	#print len(li)
	print li[7],
	if li[7]=='':
		print 'yes'

#I will have a test
# This is a test

I will have it
How do you think about it 
 I'm just editing it.
 editing
 It works very well
 heihie 
 It works well!
 I love it
 Try it again!
 It works very well