from treeFunctions import *

t = readTree( "TTGamma_V03.06_tree.root", "photonTree")

numbers = []

for e in t:
	a = e.runNumber, e.eventNumber,e.luminosityBlockNumber
	if a in numbers:
		print a, "is there several times"
	else:
		numbers.append( a)

