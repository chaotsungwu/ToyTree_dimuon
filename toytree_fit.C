#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
double binW;
//exponential background function
Double_t EXP(Double_t *x,Double_t *par)
{
	Double_t composeE1 = par[1]*x[0];
	Double_t composeE2 = par[0]*TMath::Exp(composeE1);
	return(composeE2);
}
//gaussian peak function
Double_t GAU(Double_t *x,Double_t *par)
{	
	Double_t compose1     = 0;   //initalize to 0
	if(par[2]!=0)compose1 =((x[0]-par[1])/par[2]);
	Double_t compose2     = par[0]*(binW/(par[2]*TMath::Power(TMath::Pi()*2,0.5)));
	Double_t compose3     = compose2 * TMath::Exp(-0.5*compose1*compose1);
	return(compose3);
}
//sum of the gussian and exponential
Double_t fitfunction(Double_t *x,Double_t *par){
	return(EXP(x,par)+GAU(x,&par[2]));//the EXP use the first two index of par, so GAU must start with the third 3
}
void toytree_fit()   //we always name the main function the name same as file
{
    TCanvas *c = new TCanvas("Results","Results",600,600);
    gStyle->SetOptStat(0);   //strip off the original numerical legend.
    gStyle->SetLegendBorderSize(0);   //Strip off the legend box
	TFile *file= new TFile("Results.root");           //ToyTree.root is the name of the file.
	TH1F *h_mass = (TH1F*)file -> Get("Mass Reconstruction");   //redirect to the graph of Mass Reconstruction in "Results.root"
	binW = h_mass->GetBinWidth(0);
	double par[5];
	double ndf =0.0;
	double chi_squr=0.0;
    double error_2 = 0;
    double error_3 = 0;
    double error_4 = 0;
	//("name",function,function down bound,function up bound,number of parameters)
	TF1 *fun = new TF1("fun",fitfunction,50,140,5);
    TF1 *fun1 = new TF1("fun1",GAU,50,130,3);
    TF1 *fun2 = new TF1("fun2",EXP,50,130,2);
    fun -> SetParLimits(0,0,100000);
    fun -> SetParLimits(1,-1,0);
    fun -> SetParLimits(2,0,1000000);
    fun -> SetParLimits(3,80,100);
    fun -> SetParLimits(4,0,65);
	h_mass -> Fit("fun","","",50,130);
    fun -> GetParameters(par);
    //get "Mean" and "chi_squr/ndf"
	chi_squr = fun -> GetChisquare();
	ndf = fun -> GetNDF();
	cout << endl;
	cout <<" Mean: " << par[3] << endl;
	cout <<" chi_squr/ndf: " << (double) chi_squr/ndf << endl;
    error_2 = fun -> GetParError(2);
    error_3 = fun -> GetParError(3);
    error_4 = fun -> GetParError(4);
    //Draw the pure signal in the mass distribution
    fun1 -> SetParameter(0,par[2]);
    fun1 -> SetParameter(1,par[3]);
    fun1 -> SetParameter(2,par[4]);
    fun1 -> SetLineColor(kGreen);
    fun1 -> SetFillColor(kGreen);
    fun1 -> SetFillStyle(3003);
    fun1 -> Draw("SAME");
    fun2 -> SetParameter(0,par[0]);
    fun2 -> SetParameter(1,par[1]);
    fun2 -> SetLineColor(kRed);
    fun2 -> SetLineStyle(4);
    fun2 -> Draw("SAME");
	//Add the axes
	TAxis *Xaxis = h_mass -> GetXaxis();
	TAxis *Yaxis = h_mass -> GetYaxis();
	Xaxis -> SetTitle("Mass(GeV/c^{2})");Xaxis -> SetTitleOffset(1.5);
	Yaxis -> SetTitle("Number of events");Yaxis -> SetTitleOffset(1.5);
    //Draw the legend
    TLegend *legend=new TLegend(0.60,0.65,0.85,0.85);
    //legend->SetTextFont(72);
    //legend->SetTextSize(0.04);
    legend->AddEntry(h_mass,"Data","lpe");
    legend->AddEntry(fun,"Function fit","l");
    legend->AddEntry(fun1,"Pure signal","f");
    legend->AddEntry(fun2,"Background","l");
    legend->AddEntry((TObject*)0,Form("N = %.3f #pm %.3f",par[2],error_2),"");
    legend->AddEntry((TObject*)0,Form("M_{#mu#mu} = %.3f #pm %.3f",par[3],error_3),"");
    legend->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f",par[4],error_4),"");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.3f",chi_squr/ndf),"");
    legend->AddEntry((TObject*)0,Form("BinW = %.1f",binW),"");
    legend->Draw("SAME");
	//Draw the histogram
	h_mass -> Draw("same");
	TFile *output = new TFile("Results_fit.root","recreate");
	c -> Write();
}
