#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TObject.h"
#include "THStack.h"
#include<string.h>
#include<vector>
#include<algorithm>
using std::string;
using namespace std;

void effplot(string filename, string name, string save, string xlabel, string suff = "", float xlow = 0., float xhigh = 200., int rebin = 2, float bottom = 0., float top = 1.3) {

	gROOT->ForceStyle();

	TFile *file1 = new TFile(Form("%s_pu35.root",filename.c_str()));
	TH1F *hist1_num = (TH1F*)file1->Get(Form("check/%snum%s",name.c_str(),suff.c_str()));
	TH1F *hist1 = (TH1F*)file1->Get(Form("check/%sden%s",name.c_str(),suff.c_str()));
	TH1F* hist1_den = (TH1F*)hist1->Clone("hist1_den");

	hist1_num->SetName("hist1_num");

	hist1_num->Sumw2();
	hist1_den->Sumw2();
	hist1_num->Rebin(rebin);
	hist1_den->Rebin(rebin);

	TGraphAsymmErrors *hEff1 = new TGraphAsymmErrors();
        hEff1->Divide(hist1_num,hist1_den,"cl=0.683 b(1,1) mode");

	hEff1->SetMarkerStyle(20);
	hEff1->SetMarkerColor(1);
	hEff1->SetLineColor(1);
	hEff1->SetLineWidth(2);
	hist1_den->SetLineColor(1);
	hist1_den->SetLineWidth(2);

	TFile *file2 = new TFile(Form("%s_pu140.root",filename.c_str()));
	TH1F *hist2_num = (TH1F*)file2->Get(Form("check/%snum%s",name.c_str(),suff.c_str()));
	TH1F *hist2 = (TH1F*)file2->Get(Form("check/%sden%s",name.c_str(),suff.c_str()));
	TH1F* hist2_den = (TH1F*)hist2->Clone("hist2_den");

	hist2_num->SetName("hist2_num");

	hist2_num->Sumw2();
	hist2_den->Sumw2();
	hist2_num->Rebin(rebin);
	hist2_den->Rebin(rebin);

	TGraphAsymmErrors *hEff2 = new TGraphAsymmErrors();
        hEff2->Divide(hist2_num,hist2_den,"cl=0.683 b(1,1) mode");

	hEff2->SetMarkerStyle(21);
	hEff2->SetMarkerColor(2);
	hEff2->SetLineColor(2);
	hEff2->SetLineWidth(2);
	hist2_den->SetLineColor(2);
	hist2_den->SetLineWidth(2);

	TFile *file3 = new TFile(Form("%s_pu200.root",filename.c_str()));
	TH1F *hist3_num = (TH1F*)file3->Get(Form("check/%snum%s",name.c_str(),suff.c_str()));
	TH1F *hist3 = (TH1F*)file3->Get(Form("check/%sden%s",name.c_str(),suff.c_str()));
	TH1F* hist3_den = (TH1F*)hist3->Clone("hist3_den");

	hist3_num->SetName("hist3_num");

	hist3_num->Sumw2();
	hist3_den->Sumw2();
	hist3_num->Rebin(rebin);
	hist3_den->Rebin(rebin);

	TGraphAsymmErrors *hEff3 = new TGraphAsymmErrors();
        hEff3->Divide(hist3_num,hist3_den,"cl=0.683 b(1,1) mode");

	hEff3->SetMarkerStyle(22);
	hEff3->SetMarkerColor(4);
	hEff3->SetLineColor(4);
	hEff3->SetLineWidth(2);
	hist3_den->SetLineColor(4);
	hist3_den->SetLineWidth(2);

	TCanvas *cst = new TCanvas("cst","cst",10,10,700,700);
	gStyle->SetOptStat(0);

	hEff1->SetName("hEff1");
	hEff2->SetName("hEff2");
	hEff3->SetName("hEff3");

	TMultiGraph *mgr = new TMultiGraph();
	mgr->Add(hEff1);
	mgr->Add(hEff2);
	mgr->Add(hEff3);

	hist1_num->GetYaxis()->SetTitle("Efficiency");
	hist1_num->GetXaxis()->SetTitle(xlabel.c_str());
	hist1_num->GetXaxis()->SetRangeUser(xlow,xhigh);
	hist1_num->SetMaximum(top);
	hist1_num->SetMinimum(bottom);
	hist1_num->Draw("axis");
	mgr->Draw("P");
	double hmax = hist3_den->GetMaximum();
	if (hmax==0) return;
	hist1_den->Scale(0.5/hmax);
	hist2_den->Scale(0.5/hmax);
	hist3_den->Scale(0.5/hmax);
	hist1_den->Draw("hist same");
	hist2_den->Draw("hist same");
	hist3_den->Draw("hist same");

	cst->Update();
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
	    gPad->GetUxmax(),gPad->GetUymax(),bottom,2.*hmax*top,510,"+L");
	axis->SetLineColor(kBlue);
	axis->SetTextColor(kBlue);
	axis->Draw();
	
	TLegend *leg = new TLegend(0.1,0.75,0.9,0.9);
	leg->SetNColumns(3);
	TLegendEntry *entry1 = leg->AddEntry("hEff1","PU35 - Efficiency","ep");
	TLegendEntry *entry2 = leg->AddEntry("hEff2","PU140 - Efficiency","ep");
	TLegendEntry *entry3 = leg->AddEntry("hEff3","PU200 - Efficiency","ep");
	TLegendEntry *entry4 = leg->AddEntry("hist1_den","PU35 - Distribution","l");
	TLegendEntry *entry5 = leg->AddEntry("hist2_den","PU140 - Distribution","l");
	TLegendEntry *entry6 = leg->AddEntry("hist3_den","PU200 - Distribution","l");
	leg->Draw();

	cst->SaveAs(Form("%sEff_%s.png",name.c_str(),save.c_str()));
	cst->SaveAs(Form("%sEff_%s.pdf",name.c_str(),save.c_str()));
}

void distplot(string filename, string name, string save, string xlabel, string suff = "", float xlow = 0., float xhigh = 200., int rebin = 2, float bottom = 0., float top = 80.) {

	gROOT->ForceStyle();

	TFile *file1 = new TFile(Form("%s_pu35.root",filename.c_str()));
	TH1F *hist1 = (TH1F*)file1->Get(Form("check/%s%s",name.c_str(),suff.c_str()));

	hist1->SetName("hist1");

	hist1->Sumw2();
	hist1->Rebin(rebin);

	hist1->SetLineColor(1);
	hist1->SetLineWidth(2);

	TFile *file2 = new TFile(Form("%s_pu140.root",filename.c_str()));
	TH1F *hist2 = (TH1F*)file2->Get(Form("check/%s%s",name.c_str(),suff.c_str()));

	hist2->SetName("hist2");

	hist2->Sumw2();
	hist2->Rebin(rebin);

	hist2->SetLineColor(2);
	hist2->SetLineWidth(2);

	TFile *file3 = new TFile(Form("%s_pu200.root",filename.c_str()));
	TH1F *hist3 = (TH1F*)file3->Get(Form("check/%s%s",name.c_str(),suff.c_str()));

	hist3->SetName("hist3");

	hist3->Sumw2();
	hist3->Rebin(rebin);

	hist3->SetLineColor(4);
	hist3->SetLineWidth(2);

	TCanvas *cst = new TCanvas("cst","cst",10,10,700,700);
	gStyle->SetOptStat(0);

	hist1->GetYaxis()->SetTitle("Arb. units");
	hist1->GetXaxis()->SetTitle(xlabel.c_str());
	hist1->GetXaxis()->SetRangeUser(xlow,xhigh);
	double hmax = hist3->GetMaximum();
	if (hmax==0) return;
	hist1->Scale(1/hist1->Integral());
	hist2->Scale(1/hist2->Integral());
	hist3->Scale(1/hist3->Integral());
	hist1->SetMaximum(top);
	//hist1->SetMinimum(bottom);
	hist1->Draw("hist");
	hist2->Draw("hist same");
	hist3->Draw("hist same");

	TLegend *leg = new TLegend(0.1,0.75,0.9,0.9);
	leg->SetNColumns(3);
	TLegendEntry *entry4 = leg->AddEntry("hist1","PU35 - Distribution","l");
	TLegendEntry *entry5 = leg->AddEntry("hist2","PU140 - Distribution","l");
	TLegendEntry *entry6 = leg->AddEntry("hist3","PU200 - Distribution","l");
	leg->Draw();

	cst->SetLogy();

	cst->SaveAs(Form("%sDist_%s.png",name.c_str(),save.c_str()));
	cst->SaveAs(Form("%sDist_%s.pdf",name.c_str(),save.c_str()));
}

void EffPlot()
{
        effplot("histo_ttbar", "el", "ttbar", "gen e p_{T}","",0.,200.);
        effplot("histo_ttbar", "el", "gen_ttbar", "gen e p_{T}","3",0.,200.);

        effplot("histo_ttbar", "eleta", "ttbar", "gen e #eta","",-3.,3.,4);
        effplot("histo_ttbar", "eleta", "gen_ttbar", "gen e #eta","3",-3.,3.,4);

        effplot("histo_ttbar", "mu", "ttbar", "gen #mu p_{T}","",0.,200.);
        effplot("histo_ttbar", "mu", "gen_ttbar", "gen #mu p_{T}","3",0.,200.);

        effplot("histo_ttbar", "mueta", "ttbar", "gen #mu #eta","",-3.,3.,4);
        effplot("histo_ttbar", "mueta", "gen_ttbar", "gen #mu #eta","3",-3.,3.,4);

        distplot("histo_ttbar", "jtpt", "ttbar", "jet p_{T}","",0.,1000.);
        distplot("histo_ttbar", "jteta", "ttbar", "jet #eta","",-5.,5.);
        distplot("histo_ttbar", "mtpt", "ttbar", "MET","",0.,1000.);
}
