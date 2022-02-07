#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
double binW;   //Can be discarded
/*
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
*/
void toytree_pt()   //we always name the main function the name same as file
{  
    gStyle->SetOptStat(0);   //strip off the original numerical legend.
    TFile *file= new TFile("ToyTree.root");           //ToyTree.root is the name of the file.
	TTreeReader reader("ToyTree",file);               //ToyTree is the name of tree in the file.
	TTreeReaderArray<double> pt1(reader,"mu1_pt");
	TTreeReaderArray<double> pt2(reader,"mu2_pt");
	TTreeReaderArray<double> eta1(reader,"mu1_eta");
	TTreeReaderArray<double> eta2(reader,"mu2_eta");
	TTreeReaderArray<double> phi1(reader,"mu1_phi");
	TTreeReaderArray<double> phi2(reader,"mu2_phi");
  	//TH1F declared.(TH1F is the histogram.) 
	//(filename,title,bin num,lower bound,upper bound) 
	TH1F *h_mass = new TH1F("Mass Reconstruction","Mass Reconstruction",200,50,150);
    TH1F *mass_signal = new TH1F("mass_signal","mass_signal",200,50,150);
    TH1F *pt_signal = new TH1F("pt_signal","pt_signal",50,0,50);   //The graph is going to be used in "toytree_final4.C"
    TH1F *mass_BG = new TH1F("mass_BG","mass_BG",200,50,150);
    TH1F *pt_BG = new TH1F("pt_BG","pt_BG",50,0,50);   //The graph is going to be used in "toytree_final4.C"
	h_mass->Sumw2();   //To tranfer the data of histogram into the data points
    pt_signal->Sumw2();
    pt_BG->Sumw2();
    binW = h_mass->GetBinWidth(0);   //Can be discarded
    //Fill in to show the invariant mass distribution
	double E = 0;
    double E_squr = 0;
    double p_squr = 0;
	double Mass = 0;
    double Mmuon = 0.1056583715;
    double pt = 0;
    double pt1_squr = 0;
    double pt2_squr = 0;
    double sinh1_squr = 0;
    double sinh2_squr = 0;
    double Mmuon_squr = 0;
    double px,py,pz = 0;
    double par[5];
    double ndf =0.0;
    double chi=0.0;
	while (reader.Next())
	{
        for(int i=0,n = pt1.GetSize();i<n;i++)
		{
            pt1_squr=pow(pt1[i],2);
            pt2_squr=pow(pt2[i],2);
            sinh1_squr=pow(sinh(eta1[i]),2);
            sinh2_squr=pow(sinh(eta2[i]),2);
            Mmuon_squr=pow(Mmuon,2);
            E=pow(pt1_squr*sinh1_squr+pt1_squr+Mmuon_squr,0.5)+pow(pt2_squr*sinh2_squr+pt2_squr+Mmuon_squr,0.5);   //Energy
            E_squr=pow(E,2);   //Energy square
            px=pt1[i]*cos(phi1[i])+pt2[i]*cos(phi2[i]);
            py=pt1[i]*sin(phi1[i])+pt2[i]*sin(phi2[i]);
            pz=pt1[i]*sinh(eta1[i])+pt2[i]*sinh(eta2[i]);
            pt=pow(pow(px,2)+pow(py,2),0.5);
            p_squr=pow(px,2)+pow(py,2)+pow(pz,2);   //Momentum square
            Mass=pow(E_squr-p_squr,0.5);   //Mass in each event
			h_mass -> Fill(Mass);   //Fill in
            if (Mass>(87.6486-3*5.12029) && Mass<(87.6486+3*5.12029)){
                mass_signal -> Fill(Mass);
                pt_signal -> Fill(pt);   //Get to do normalization
            }
            if ((50<Mass && Mass<(87.6486-3*5.12029)) || (130>Mass && Mass>(87.6486+3*5.12029))){   //Because of the lower bound 50 and upper bound 130
                mass_BG -> Fill(Mass);
                pt_BG -> Fill(pt);   //Get to do normalization
            }
		}
	}
/*
    //("name",function,function down bound,function up bound,number of parameters)
    TF1 *fun = new TF1("fun",fitfunction,50,140,5);
    fun -> SetParLimits(0,0,100000);
    fun -> SetParLimits(1,-1,0);
    fun -> SetParLimits(2,0,1000000);
    fun -> SetParLimits(3,80,100);
    fun -> SetParLimits(4,0,100000);
    h_mass -> Fit("fun","","",50,130);
    fun -> GetParameters(par);
    //get "Mean" and "chi_squr/ndf"
    chi = fun -> GetChisquare();
    ndf = fun -> GetNDF();
    cout << endl;
    cout <<" Mean: " << par[3] << endl;
    cout <<" chi_squr/ndf: " << (double) chi/ndf << endl;
    //Add the axes
    TAxis *Xaxis = h_mass -> GetXaxis();
    TAxis *Yaxis = h_mass -> GetYaxis();
    Xaxis -> SetTitle("Mass(GeV)");Xaxis -> SetTitleOffset(1.5);
    Yaxis -> SetTitle("Number of events");Yaxis -> SetTitleOffset(1.5);
*/
    //Draw the histogram
	TCanvas *a = new TCanvas("Results","Results",600,600);
	h_mass -> Draw();
    TCanvas *b = new TCanvas("mass_signal","mass_signal",600,600);
    mass_signal -> Draw();
    TCanvas *c = new TCanvas("pt_signal","pt_signal",600,600);
    pt_signal -> Draw();
    TCanvas *d = new TCanvas("mass_BG","mass_BG",600,600);
    mass_BG -> Draw();
    TCanvas *e = new TCanvas("pt_BG","pt_BG",600,600);
    pt_BG -> Draw();
    //Write the graphs into "Results.root"
    TFile *output= new TFile("Results.root","recreate");
    h_mass -> Write();
    mass_signal -> Write();
    pt_signal -> Write();
    mass_BG -> Write();
    pt_BG -> Write();
    //Draw the legend
    //TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
    //legend->SetTextFont(72);
    //legend->SetTextSize(0.04);
    //legend->AddEntry(h_mass,"Data","lpe");
    //legend->AddEntry(fun,"Function fit","l");
    //legend->Draw("same");
}
