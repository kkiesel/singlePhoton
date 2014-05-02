
from treeFunctions import *

"""
filenames = ["B_1300_1720_375_V03.07_tree.root", "B_1700_1120_375_V03.07_tree.root", "W_1700_720_375_V03.07_tree.root", "W_900_1720_375_V03.07_tree.root"]
acc = {"B_1300_1720_375_V03.07_tree.root": (7324.1353999, 4336.9018),
  "B_1700_1120_375_V03.07_tree.root":(6557.3294199999, 2473.5701),
  "W_1700_720_375_V03.07_tree.root":( 4409.4644,847.6713),
  "W_900_1720_375_V03.07_tree.root":(4340.91556, 1555.309399999)}

for f in filenames:
	n = 100 if f[0] == "B" else 600


	t = readTree( f, "photonTree" )
	withoutCutInclusive = 1.*t.Draw("met", "met>100", "goff" )/n
	withoutCut          = 1.*t.Draw("met", "met>350", "goff")/n
	htCutInclusive = 1.*t.Draw("met", "met>100 && ht>500", "goff")/n
	htCut          = 1.*t.Draw("met", "met>350 && ht>500", "goff")/n
	ptCutInclusive = 1.*t.Draw("met", "met>100 && photons[0].ptJet()>110", "goff")/n
	ptCut          = 1.*t.Draw("met", "met>350 && photons[0].ptJet()>110", "goff")/n
	leptonCutInclusive = 1.*t.Draw("met", "met>100 && !@electrons.size() && !@muons.size()", "goff")/n
	leptonCut          = 1.*t.Draw("met", "met>350 && !@electrons.size() && !@muons.size()", "goff")/n
	allCutInclusive = 1.*t.Draw("met", "met>100 && photons[0].ptJet()>110 && ht>500 && !@electrons.size() && !@muons.size()", "goff")/n
	allCut          = 1.*t.Draw("met", "met>350 && photons[0].ptJet()>110 && ht>500 && !@electrons.size() && !@muons.size()", "goff")/n

	oldInclusive, old = acc[f]

	oldInclusive /= n
	old /= n

	print "\n",f
	print (old-withoutCut)/withoutCut*100
	print (withoutCut-withoutCut)/withoutCut*100
	print (htCut-withoutCut)/withoutCut*100
	print (ptCut-withoutCut)/withoutCut*100
	print (leptonCut-withoutCut)/withoutCut*100
	print (allCut-withoutCut)/withoutCut*100
"""



	#print "all cuts exc. lepton veto ",
	#print 1.*t.Draw("met", "met>100 && photons[0].ptJet()>110 && ht>500", "goff")/n, 1.*t.Draw("met", "met>350 && photons[0].ptJet()>110 && ht>500", "goff")/n


# calculate significance
# use last met bin met>350

def significanceAsimov( s, b ):
	import math
	return math.sqrt( 2* ( (s+b)*math.log( 1+s/b ) - s ) )

def significanceBkg( s, b, s_b ):
	import math
	return s/math.sqrt( b + s_b**2 )

if __name__ == "__main__":
# signal acceptances:
	b1300 = { "nCut": 42, "htCut": 42, "ptCut": 40, "leptonCut": 34, "allCut": 32, "xSec":5.33e-3 }
	b1700 = { "nCut": 25, "htCut": 24, "ptCut": 23, "leptonCut": 15, "allCut": 14, "xSec":9.7e-3  }
	w1700 = { "nCut": 1.3, "htCut": 1.3, "ptCut": 1.2, "leptonCut": 0.75, "allCut": 0.66, "xSec":0.316 }
	w_900 = { "nCut": 2.4, "htCut": 2.4, "ptCut": 2.2, "leptonCut": 1.7, "allCut": 1.6, "xSec":8.37e-2 }

	signals = { "b1300": b1300, "b1700": b1700, "w1700": w1700, "w_900": w_900 }

#selected events
	qcd = { "nCut": 6.288, "htCut": 6.0376, "ptCut": 0, "leptonCut": 4.4478, "allCut": 4.27 }

	lumi = 3932

	for cut, sel in qcd.iteritems():
		print
		print cut
		for spoint, signal in signals.iteritems():
			print spoint,
			s = signal[cut]*signal["xSec"]*lumi/100
			print significanceAsimov( s, sel ) if sel else "Warning: 0 background"

