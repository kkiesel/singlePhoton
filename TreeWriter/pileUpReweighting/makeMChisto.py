import ROOT

if __name__ == "__main__":

	histo = ROOT.TH1F("pileup","pileup",60,0,60)
	Summer2012_S10 = [
				2.560E-06,
				5.239E-06,
				1.420E-05,
				5.005E-05,
				1.001E-04,
				2.705E-04,
				1.999E-03,
				6.097E-03,
				1.046E-02,
				1.383E-02,
				1.685E-02,
				2.055E-02,
				2.572E-02,
				3.262E-02,
				4.121E-02,
				4.977E-02,
				5.539E-02,
				5.725E-02,
				5.607E-02,
				5.312E-02,
				5.008E-02,
				4.763E-02,
				4.558E-02,
				4.363E-02,
				4.159E-02,
				3.933E-02,
				3.681E-02,
				3.406E-02,
				3.116E-02,
				2.818E-02,
				2.519E-02,
				2.226E-02,
				1.946E-02,
				1.682E-02,
				1.437E-02,
				1.215E-02,
				1.016E-02,
				8.400E-03,
				6.873E-03,
				5.564E-03,
				4.457E-03,
				3.533E-03,
				2.772E-03,
				2.154E-03,
				1.656E-03,
				1.261E-03,
				9.513E-04,
				7.107E-04,
				5.259E-04,
				3.856E-04,
				2.801E-04,
				2.017E-04,
				1.439E-04,
				1.017E-04,
				7.126E-05,
				4.948E-05,
				3.405E-05,
				2.322E-05,
				1.570E-05,
				5.005E-06]

	for i in range(1,61):
		histo.SetBinContent(i,Summer2012_S10[i-1]);


	datafile = ROOT.TFile("PU_dist.root")
	dataHist = datafile.Get("pileup")

	for h in [dataHist, histo]:
		h.Scale( 1./h.Integral() )

	dataHist.Divide( histo )

	weightFile = ROOT.TFile("puWeights.root", "recreate")
	weightFile.cd()
	dataHist.Write()
	weightFile.Close()

	file = ROOT.TFile("mc_PU_dist_S10.root","recreate")
	file.cd()
	histo.Write()
	file.Close()

	datafile.Close()
