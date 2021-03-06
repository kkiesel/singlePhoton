{
//=========Macro generated from canvas: limit_plot/limit
//=========  (Wed Jul 16 14:32:09 2014) by ROOT version5.34/03
   TCanvas *limit_plot = new TCanvas("limit_plot", "limit",1617,44,1000,900);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   limit_plot->SetHighLightColor(2);
   limit_plot->Range(0,0,1,1);
   limit_plot->SetFillColor(0);
   limit_plot->SetBorderMode(0);
   limit_plot->SetBorderSize(2);
   limit_plot->SetTickx(1);
   limit_plot->SetTicky(1);
   limit_plot->SetLeftMargin(0.15);
   limit_plot->SetBottomMargin(0.2);
   limit_plot->SetFrameFillStyle(0);
   limit_plot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: top_pad
   TPad *top_pad = new TPad("top_pad", "top_pad",0,0.8,1,1);
   top_pad->Draw();
   top_pad->cd();
   top_pad->Range(0,0,1,1);
   top_pad->SetFillColor(0);
   top_pad->SetBorderMode(0);
   top_pad->SetBorderSize(2);
   top_pad->SetTickx(1);
   top_pad->SetTicky(1);
   top_pad->SetLeftMargin(0.15);
   top_pad->SetTopMargin(0);
   top_pad->SetBottomMargin(0.4);
   top_pad->SetFrameFillStyle(0);
   top_pad->SetFrameBorderMode(0);
   
   TPaveText *pt = new TPaveText(0.1495984,0.01376147,0.9006024,0.9541284,"brNDC");
   pt->SetFillColor(0);
   pt->SetLineWidth(3);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.172);
   TText *text = pt->AddText("pp #rightarrow #tilde{g}#tilde{g},    #tilde{g} #rightarrow q#bar{q} + #tilde{#chi}^{0}_{1},    #tilde{#chi}^{0}_{1} #rightarrow #tilde{G} #gamma");

   text = pt->AddText("CMS Simulation #sqrt{s} = 8 TeV");
   pt->Draw();
   top_pad->Modified();
   limit_plot->cd();
  
// ------------>Primitives in pad: low_bad
   low_bad = new TPad("low_bad", "low_pad",0,0,1,0.8);
   low_bad->Draw();
   low_bad->cd();
   low_bad->Range(-0.2871422,-3.254342,5.643016,2.750297);
   low_bad->SetFillColor(0);
   low_bad->SetBorderMode(0);
   low_bad->SetBorderSize(2);
   low_bad->SetLogy();
   low_bad->SetTickx(1);
   low_bad->SetTicky(1);
   low_bad->SetLeftMargin(0.1495984);
   low_bad->SetTopMargin(0);
   low_bad->SetBottomMargin(0.2006881);
   low_bad->SetFrameFillStyle(0);
   low_bad->SetFrameBorderMode(0);
   low_bad->SetFrameFillStyle(0);
   low_bad->SetFrameBorderMode(0);
   Double_t xAxis1[90] = {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.05}; 
   
   TH1F *hist_high = new TH1F("hist_high","High Rsq",89, xAxis1);
   hist_high->SetBinContent(1,498.1906);
   hist_high->SetBinContent(2,315.1296);
   hist_high->SetBinContent(3,205.329);
   hist_high->SetBinContent(4,137.3048);
   hist_high->SetBinContent(5,93.94713);
   hist_high->SetBinContent(6,65.60689);
   hist_high->SetBinContent(7,46.66177);
   hist_high->SetBinContent(8,33.73918);
   hist_high->SetBinContent(9,24.76259);
   hist_high->SetBinContent(10,18.42317);
   hist_high->SetBinContent(11,13.87818);
   hist_high->SetBinContent(12,10.57444);
   hist_high->SetBinContent(13,8.142314);
   hist_high->SetBinContent(14,6.330793);
   hist_high->SetBinContent(15,4.96684);
   hist_high->SetBinContent(16,3.929515);
   hist_high->SetBinContent(17,3.133195);
   hist_high->SetBinContent(18,2.516536);
   hist_high->SetBinContent(19,2.035096);
   hist_high->SetBinContent(20,1.656344);
   hist_high->SetBinContent(21,1.356231);
   hist_high->SetBinContent(22,1.116819);
   hist_high->SetBinContent(23,0.9246118);
   hist_high->SetBinContent(24,0.769372);
   hist_high->SetBinContent(25,0.6432748);
   hist_high->SetBinContent(26,0.5402965);
   hist_high->SetBinContent(27,0.4557676);
   hist_high->SetBinContent(28,0.3860454);
   hist_high->SetBinContent(29,0.32827);
   hist_high->SetBinContent(30,0.2801834);
   hist_high->SetBinContent(31,0.2399928);
   hist_high->SetBinContent(32,0.2062667);
   hist_high->SetBinContent(33,0.1778571);
   hist_high->SetBinContent(34,0.153838);
   hist_high->SetBinContent(35,0.1334598);
   hist_high->SetBinContent(36,0.1161123);
   hist_high->SetBinContent(37,0.1012971);
   hist_high->SetBinContent(38,0.08860548);
   hist_high->SetBinContent(39,0.07770059);
   hist_high->SetBinContent(40,0.0683041);
   hist_high->SetBinContent(41,0.06018505);
   hist_high->SetBinContent(42,0.05315118);
   hist_high->SetBinContent(43,0.04704184);
   hist_high->SetBinContent(44,0.04172244);
   hist_high->SetBinContent(45,0.0370798);
   hist_high->SetBinContent(46,0.03301853);
   hist_high->SetBinContent(47,0.02945793);
   hist_high->SetBinContent(48,0.02632957);
   hist_high->SetBinContent(49,0.02357528);
   hist_high->SetBinContent(50,0.02114544);
   hist_high->SetBinContent(51,0.01899767);
   hist_high->SetBinContent(52,0.01709564);
   hist_high->SetBinContent(53,0.01540816);
   hist_high->SetBinContent(54,0.01390837);
   hist_high->SetBinContent(55,0.01257309);
   hist_high->SetBinContent(56,0.0113823);
   hist_high->SetBinContent(57,0.01031864);
   hist_high->SetBinContent(58,0.009367035);
   hist_high->SetBinContent(59,0.00851438);
   hist_high->SetBinContent(60,0.007749247);
   hist_high->SetBinContent(61,0.007061657);
   hist_high->SetBinContent(62,0.006442877);
   hist_high->SetBinContent(63,0.00588526);
   hist_high->SetBinContent(64,0.005382085);
   hist_high->SetBinContent(65,0.004927447);
   hist_high->SetBinContent(66,0.004516141);
   hist_high->SetBinContent(67,0.004143578);
   hist_high->SetBinContent(68,0.0038057);
   hist_high->SetBinContent(69,0.003498917);
   hist_high->SetBinContent(70,0.003220048);
   hist_high->SetBinContent(71,0.002966268);
   hist_high->SetBinContent(72,0.002735067);
   hist_high->SetBinContent(73,0.002524212);
   hist_high->SetBinContent(74,0.002331709);
   hist_high->SetBinContent(75,0.002155784);
   hist_high->SetBinContent(76,0.001994847);
   hist_high->SetBinContent(77,0.001847478);
   hist_high->SetBinContent(78,0.001712405);
   hist_high->SetBinContent(79,0.001588485);
   hist_high->SetBinContent(80,0.001474694);
   hist_high->SetBinContent(81,0.001370111);
   hist_high->SetBinContent(82,0.001273907);
   hist_high->SetBinContent(83,0.001185334);
   hist_high->SetBinContent(84,0.001103718);
   hist_high->SetBinContent(85,0.001028452);
   hist_high->SetBinContent(86,0.0009589846);
   hist_high->SetBinContent(87,0.0008948195);
   hist_high->SetBinContent(88,0.0008355058);
   hist_high->SetBinContent(89,0.0007806351);
   hist_high->SetMinimum(0.008927256);
   hist_high->SetMaximum(562.7262);
   hist_high->SetEntries(89);
   hist_high->SetStats(0);
   hist_high->SetLineColor(13);
   hist_high->SetLineWidth(3);
   hist_high->GetXaxis()->SetTitle("M_{R} (TeV)");
   hist_high->GetXaxis()->SetLabelFont(42);
   hist_high->GetXaxis()->SetLabelSize(0.05);
   hist_high->GetXaxis()->SetTitleSize(0.065);
   hist_high->GetXaxis()->SetTitleFont(42);
   hist_high->GetYaxis()->SetTitle("N Events");
   hist_high->GetYaxis()->SetLabelFont(42);
   hist_high->GetYaxis()->SetLabelSize(0.05);
   hist_high->GetYaxis()->SetTitleSize(0.065);
   hist_high->GetYaxis()->SetTitleFont(42);
   hist_high->GetZaxis()->SetLabelFont(42);
   hist_high->GetZaxis()->SetLabelSize(0.035);
   hist_high->GetZaxis()->SetTitleSize(0.035);
   hist_high->GetZaxis()->SetTitleFont(42);
   hist_high->Draw("");
   Double_t xAxis2[90] = {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.05}; 
   
   TH1F *hist_high_sig1 = new TH1F("hist_high_sig1","hist_high_sig1",89, xAxis2);
   hist_high_sig1->SetBinContent(1,0.012795);
   hist_high_sig1->SetBinContent(2,0.02559);
   hist_high_sig1->SetBinContent(3,0.012795);
   hist_high_sig1->SetBinContent(4,0.02559);
   hist_high_sig1->SetBinContent(5,0.015354);
   hist_high_sig1->SetBinContent(6,0.020472);
   hist_high_sig1->SetBinContent(7,0.05885699);
   hist_high_sig1->SetBinContent(8,0.07932898);
   hist_high_sig1->SetBinContent(9,0.06909299);
   hist_high_sig1->SetBinContent(10,0.09468298);
   hist_high_sig1->SetBinContent(11,0.122832);
   hist_high_sig1->SetBinContent(12,0.1714531);
   hist_high_sig1->SetBinContent(13,0.2456643);
   hist_high_sig1->SetBinContent(14,0.2661363);
   hist_high_sig1->SetBinContent(15,0.2942854);
   hist_high_sig1->SetBinContent(16,0.3352295);
   hist_high_sig1->SetBinContent(17,0.3429065);
   hist_high_sig1->SetBinContent(18,0.4043226);
   hist_high_sig1->SetBinContent(19,0.5322729);
   hist_high_sig1->SetBinContent(20,0.4990059);
   hist_high_sig1->SetBinContent(21,0.5092419);
   hist_high_sig1->SetBinContent(22,0.547627);
   hist_high_sig1->SetBinContent(23,0.5860121);
   hist_high_sig1->SetBinContent(24,0.547627);
   hist_high_sig1->SetBinContent(25,0.578335);
   hist_high_sig1->SetBinContent(26,0.6192791);
   hist_high_sig1->SetBinContent(27,0.5015649);
   hist_high_sig1->SetBinContent(28,0.4836518);
   hist_high_sig1->SetBinContent(29,0.4580618);
   hist_high_sig1->SetBinContent(30,0.4503847);
   hist_high_sig1->SetBinContent(31,0.4068816);
   hist_high_sig1->SetBinContent(32,0.3684965);
   hist_high_sig1->SetBinContent(33,0.2994034);
   hist_high_sig1->SetBinContent(34,0.2431052);
   hist_high_sig1->SetBinContent(35,0.1868071);
   hist_high_sig1->SetBinContent(36,0.143304);
   hist_high_sig1->SetBinContent(37,0.1714531);
   hist_high_sig1->SetBinContent(38,0.138186);
   hist_high_sig1->SetBinContent(39,0.110037);
   hist_high_sig1->SetBinContent(40,0.112596);
   hist_high_sig1->SetBinContent(41,0.07165199);
   hist_high_sig1->SetBinContent(42,0.08700598);
   hist_high_sig1->SetBinContent(43,0.05629799);
   hist_high_sig1->SetBinContent(44,0.05629799);
   hist_high_sig1->SetBinContent(45,0.028149);
   hist_high_sig1->SetBinContent(46,0.035826);
   hist_high_sig1->SetBinContent(47,0.020472);
   hist_high_sig1->SetBinContent(48,0.028149);
   hist_high_sig1->SetBinContent(49,0.028149);
   hist_high_sig1->SetBinContent(50,0.007677);
   hist_high_sig1->SetBinContent(51,0.010236);
   hist_high_sig1->SetBinContent(52,0.007677);
   hist_high_sig1->SetBinContent(53,0.015354);
   hist_high_sig1->SetBinContent(54,0.010236);
   hist_high_sig1->SetBinContent(55,0.007677);
   hist_high_sig1->SetBinContent(56,0.007677);
   hist_high_sig1->SetBinContent(58,0.002559);
   hist_high_sig1->SetBinContent(59,0.002559);
   hist_high_sig1->SetBinContent(61,0.002559);
   hist_high_sig1->SetBinContent(62,0.005118);
   hist_high_sig1->SetBinContent(64,0.002559);
   hist_high_sig1->SetMinimum(0.008927256);
   hist_high_sig1->SetMaximum(562.7262);
   hist_high_sig1->SetEntries(4527);
   hist_high_sig1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff9999");
   hist_high_sig1->SetLineColor(ci);
   hist_high_sig1->SetLineWidth(3);
   hist_high_sig1->GetXaxis()->SetLabelFont(42);
   hist_high_sig1->GetXaxis()->SetLabelSize(0.035);
   hist_high_sig1->GetXaxis()->SetTitleSize(0.035);
   hist_high_sig1->GetXaxis()->SetTitleFont(42);
   hist_high_sig1->GetYaxis()->SetLabelFont(42);
   hist_high_sig1->GetYaxis()->SetLabelSize(0.035);
   hist_high_sig1->GetYaxis()->SetTitleSize(0.035);
   hist_high_sig1->GetYaxis()->SetTitleFont(42);
   hist_high_sig1->GetZaxis()->SetLabelFont(42);
   hist_high_sig1->GetZaxis()->SetLabelSize(0.035);
   hist_high_sig1->GetZaxis()->SetTitleSize(0.035);
   hist_high_sig1->GetZaxis()->SetTitleFont(42);
   hist_high_sig1->Draw("same");
   Double_t xAxis3[90] = {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.05}; 
   
   TH1F *hist_high_sig2 = new TH1F("hist_high_sig2","hist_high_sig2",89, xAxis3);
   hist_high_sig2->SetBinContent(1,0.043503);
   hist_high_sig2->SetBinContent(2,0.05117999);
   hist_high_sig2->SetBinContent(3,0.07676999);
   hist_high_sig2->SetBinContent(4,0.09468298);
   hist_high_sig2->SetBinContent(5,0.112596);
   hist_high_sig2->SetBinContent(6,0.138186);
   hist_high_sig2->SetBinContent(7,0.148422);
   hist_high_sig2->SetBinContent(8,0.2200742);
   hist_high_sig2->SetBinContent(9,0.2303102);
   hist_high_sig2->SetBinContent(10,0.3045214);
   hist_high_sig2->SetBinContent(11,0.3224344);
   hist_high_sig2->SetBinContent(12,0.3966456);
   hist_high_sig2->SetBinContent(13,0.4759748);
   hist_high_sig2->SetBinContent(14,0.5348319);
   hist_high_sig2->SetBinContent(15,0.560422);
   hist_high_sig2->SetBinContent(16,0.53995);
   hist_high_sig2->SetBinContent(17,0.6448692);
   hist_high_sig2->SetBinContent(18,0.5988071);
   hist_high_sig2->SetBinContent(19,0.575776);
   hist_high_sig2->SetBinContent(20,0.6295152);
   hist_high_sig2->SetBinContent(21,0.6243972);
   hist_high_sig2->SetBinContent(22,0.6320742);
   hist_high_sig2->SetBinContent(23,0.6141611);
   hist_high_sig2->SetBinContent(24,0.5066829);
   hist_high_sig2->SetBinContent(25,0.56554);
   hist_high_sig2->SetBinContent(26,0.4938878);
   hist_high_sig2->SetBinContent(27,0.4913288);
   hist_high_sig2->SetBinContent(28,0.4017636);
   hist_high_sig2->SetBinContent(29,0.3761736);
   hist_high_sig2->SetBinContent(30,0.4299127);
   hist_high_sig2->SetBinContent(31,0.3915276);
   hist_high_sig2->SetBinContent(32,0.2533413);
   hist_high_sig2->SetBinContent(33,0.2533413);
   hist_high_sig2->SetBinContent(34,0.2584593);
   hist_high_sig2->SetBinContent(35,0.2226332);
   hist_high_sig2->SetBinContent(36,0.2763723);
   hist_high_sig2->SetBinContent(37,0.1740121);
   hist_high_sig2->SetBinContent(38,0.138186);
   hist_high_sig2->SetBinContent(39,0.158658);
   hist_high_sig2->SetBinContent(40,0.130509);
   hist_high_sig2->SetBinContent(41,0.120273);
   hist_high_sig2->SetBinContent(42,0.08956498);
   hist_high_sig2->SetBinContent(43,0.09724198);
   hist_high_sig2->SetBinContent(44,0.08956498);
   hist_high_sig2->SetBinContent(45,0.043503);
   hist_high_sig2->SetBinContent(46,0.06141599);
   hist_high_sig2->SetBinContent(47,0.05885699);
   hist_high_sig2->SetBinContent(48,0.040944);
   hist_high_sig2->SetBinContent(49,0.035826);
   hist_high_sig2->SetBinContent(50,0.028149);
   hist_high_sig2->SetBinContent(51,0.020472);
   hist_high_sig2->SetBinContent(52,0.023031);
   hist_high_sig2->SetBinContent(53,0.015354);
   hist_high_sig2->SetBinContent(54,0.028149);
   hist_high_sig2->SetBinContent(55,0.015354);
   hist_high_sig2->SetBinContent(56,0.020472);
   hist_high_sig2->SetBinContent(57,0.005118);
   hist_high_sig2->SetBinContent(58,0.015354);
   hist_high_sig2->SetBinContent(59,0.012795);
   hist_high_sig2->SetBinContent(60,0.010236);
   hist_high_sig2->SetBinContent(61,0.015354);
   hist_high_sig2->SetBinContent(62,0.007677);
   hist_high_sig2->SetBinContent(63,0.005118);
   hist_high_sig2->SetBinContent(64,0.002559);
   hist_high_sig2->SetBinContent(66,0.007677);
   hist_high_sig2->SetBinContent(67,0.005118);
   hist_high_sig2->SetBinContent(69,0.002559);
   hist_high_sig2->SetBinContent(74,0.002559);
   hist_high_sig2->SetBinContent(83,0.002559);
   hist_high_sig2->SetMinimum(0.008927256);
   hist_high_sig2->SetMaximum(562.7262);
   hist_high_sig2->SetEntries(5852);
   hist_high_sig2->SetStats(0);

   ci = TColor::GetColor("#99ff99");
   hist_high_sig2->SetLineColor(ci);
   hist_high_sig2->SetLineWidth(3);
   hist_high_sig2->GetXaxis()->SetLabelFont(42);
   hist_high_sig2->GetXaxis()->SetLabelSize(0.035);
   hist_high_sig2->GetXaxis()->SetTitleSize(0.035);
   hist_high_sig2->GetXaxis()->SetTitleFont(42);
   hist_high_sig2->GetYaxis()->SetLabelFont(42);
   hist_high_sig2->GetYaxis()->SetLabelSize(0.035);
   hist_high_sig2->GetYaxis()->SetTitleSize(0.035);
   hist_high_sig2->GetYaxis()->SetTitleFont(42);
   hist_high_sig2->GetZaxis()->SetLabelFont(42);
   hist_high_sig2->GetZaxis()->SetLabelSize(0.035);
   hist_high_sig2->GetZaxis()->SetTitleSize(0.035);
   hist_high_sig2->GetZaxis()->SetTitleFont(42);
   hist_high_sig2->Draw("same");
   Double_t xAxis4[90] = {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.05}; 
   
   TH1F *hist_high_sig3 = new TH1F("hist_high_sig3","hist_high_sig3",89, xAxis4);
   hist_high_sig3->SetBinContent(1,0.010236);
   hist_high_sig3->SetBinContent(2,0.002559);
   hist_high_sig3->SetBinContent(3,0.007677);
   hist_high_sig3->SetBinContent(4,0.033267);
   hist_high_sig3->SetBinContent(5,0.043503);
   hist_high_sig3->SetBinContent(6,0.05373899);
   hist_high_sig3->SetBinContent(7,0.043503);
   hist_high_sig3->SetBinContent(8,0.10236);
   hist_high_sig3->SetBinContent(9,0.133068);
   hist_high_sig3->SetBinContent(10,0.1663351);
   hist_high_sig3->SetBinContent(11,0.2098382);
   hist_high_sig3->SetBinContent(12,0.2712543);
   hist_high_sig3->SetBinContent(13,0.2559003);
   hist_high_sig3->SetBinContent(14,0.3838506);
   hist_high_sig3->SetBinContent(15,0.4119996);
   hist_high_sig3->SetBinContent(16,0.5118009);
   hist_high_sig3->SetBinContent(17,0.6269562);
   hist_high_sig3->SetBinContent(18,0.5885711);
   hist_high_sig3->SetBinContent(19,0.6320742);
   hist_high_sig3->SetBinContent(20,0.6858133);
   hist_high_sig3->SetBinContent(21,0.7574655);
   hist_high_sig3->SetBinContent(22,0.6423102);
   hist_high_sig3->SetBinContent(23,0.6883723);
   hist_high_sig3->SetBinContent(24,0.6090431);
   hist_high_sig3->SetBinContent(25,0.575776);
   hist_high_sig3->SetBinContent(26,0.5860121);
   hist_high_sig3->SetBinContent(27,0.4990059);
   hist_high_sig3->SetBinContent(28,0.4938878);
   hist_high_sig3->SetBinContent(29,0.4555027);
   hist_high_sig3->SetBinContent(30,0.4196767);
   hist_high_sig3->SetBinContent(31,0.3531425);
   hist_high_sig3->SetBinContent(32,0.2763723);
   hist_high_sig3->SetBinContent(33,0.2354282);
   hist_high_sig3->SetBinContent(34,0.1996021);
   hist_high_sig3->SetBinContent(35,0.1970431);
   hist_high_sig3->SetBinContent(36,0.1893661);
   hist_high_sig3->SetBinContent(37,0.138186);
   hist_high_sig3->SetBinContent(38,0.117714);
   hist_high_sig3->SetBinContent(39,0.115155);
   hist_high_sig3->SetBinContent(40,0.08444698);
   hist_high_sig3->SetBinContent(41,0.08188798);
   hist_high_sig3->SetBinContent(42,0.048621);
   hist_high_sig3->SetBinContent(43,0.05117999);
   hist_high_sig3->SetBinContent(44,0.048621);
   hist_high_sig3->SetBinContent(45,0.046062);
   hist_high_sig3->SetBinContent(46,0.040944);
   hist_high_sig3->SetBinContent(47,0.035826);
   hist_high_sig3->SetBinContent(48,0.02559);
   hist_high_sig3->SetBinContent(49,0.023031);
   hist_high_sig3->SetBinContent(50,0.017913);
   hist_high_sig3->SetBinContent(51,0.012795);
   hist_high_sig3->SetBinContent(52,0.020472);
   hist_high_sig3->SetBinContent(53,0.010236);
   hist_high_sig3->SetBinContent(54,0.010236);
   hist_high_sig3->SetBinContent(55,0.010236);
   hist_high_sig3->SetBinContent(56,0.002559);
   hist_high_sig3->SetBinContent(57,0.005118);
   hist_high_sig3->SetBinContent(58,0.007677);
   hist_high_sig3->SetBinContent(59,0.002559);
   hist_high_sig3->SetBinContent(60,0.005118);
   hist_high_sig3->SetBinContent(61,0.002559);
   hist_high_sig3->SetBinContent(67,0.002559);
   hist_high_sig3->SetBinContent(68,0.002559);
   hist_high_sig3->SetBinContent(73,0.002559);
   hist_high_sig3->SetMinimum(0.008927256);
   hist_high_sig3->SetMaximum(562.7262);
   hist_high_sig3->SetEntries(5207);
   hist_high_sig3->SetStats(0);

   ci = TColor::GetColor("#9999ff");
   hist_high_sig3->SetLineColor(ci);
   hist_high_sig3->SetLineWidth(3);
   hist_high_sig3->GetXaxis()->SetLabelFont(42);
   hist_high_sig3->GetXaxis()->SetLabelSize(0.035);
   hist_high_sig3->GetXaxis()->SetTitleSize(0.035);
   hist_high_sig3->GetXaxis()->SetTitleFont(42);
   hist_high_sig3->GetYaxis()->SetLabelFont(42);
   hist_high_sig3->GetYaxis()->SetLabelSize(0.035);
   hist_high_sig3->GetYaxis()->SetTitleSize(0.035);
   hist_high_sig3->GetYaxis()->SetTitleFont(42);
   hist_high_sig3->GetZaxis()->SetLabelFont(42);
   hist_high_sig3->GetZaxis()->SetLabelSize(0.035);
   hist_high_sig3->GetZaxis()->SetTitleSize(0.035);
   hist_high_sig3->GetZaxis()->SetTitleFont(42);
   hist_high_sig3->Draw("same");
   
   TLegend *leg = new TLegend(0.3664659,0.5174885,0.876506,0.9690367,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextSize(0.03);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","Razor #gamma#gamma + #geq 1 jet","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hist_high_sig1","m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} = 225 GeV","l");

   ci = TColor::GetColor("#ff9999");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hist_high_sig3","m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} = 675 GeV","l");

   ci = TColor::GetColor("#9999ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hist_high_sig2","m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} = 1275 GeV","l");

   ci = TColor::GetColor("#99ff99");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hist_high","Background Model","l");
   entry->SetLineColor(13);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   low_bad->Modified();
   limit_plot->cd();
   limit_plot->Modified();
   limit_plot->cd();
   limit_plot->SetSelected(limit_plot);
}
