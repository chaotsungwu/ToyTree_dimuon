#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
double binW;
double factor = 0;
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
void toytree_final()   //we always name the main function the name same as file
{  
	gStyle->SetOptStat();   //strip off the original numerical legend.
    gStyle->SetLegendBorderSize(0);   //Strip off the legend box
	TFile *file= new TFile("Results.root");           //Results.root is the name of the file.
	TH1F *pt_signal = (TH1F*)file -> Get("pt_signal");   //redirect to the graph in "Results.root"
	TH1F *pt_BG = (TH1F*)file -> Get("pt_BG");
    //pt_signal -> Sumw2();   //Can be discarded
    //pt_BG -> Sumw2();   //Can be discarded
	binW = pt_signal->GetBinWidth(0);
    double A = 0;
    double B = 0;
    double C = 0;
    double l,u = 0;
    double area1 = 0;
    double area2 = 0;
    l=87.6488-3*5.12017;
    u=87.6488+3*5.12017;
	//("name",function,function down bound,function up bound,number of parameters)
	TF1 *fun1 = new TF1("Exponential",EXP,50,130,2);
    TF1 *fun2 = new TF1("Gaussian",GAU,50,130,3);
    TCanvas *c = new TCanvas("function","function",600,600);
    fun2 -> Draw();
    fun2 -> SetLineColor(kBlue);
    fun1 -> Draw("SAME");
    fun2 -> SetTitle("Functions");
    TLegend *legend1=new TLegend(0.6,0.65,0.88,0.85);
    legend1->AddEntry(fun1,"Exponential","l");
    legend1->AddEntry(fun2,"Gaussian","l");
    legend1->Draw("same");
	fun1 -> SetParameter(0,15411.2);
	fun1 -> SetParameter(1,-0.0398975);
	fun2 -> SetParameter(0,374827);
	fun2 -> SetParameter(1,87.6488);
	fun2 -> SetParameter(2,5.12017);
    fun1 -> FixParameter(0,15411.2);
    fun1 -> FixParameter(1,-0.0398975);
    fun2 -> FixParameter(0,374827);
    fun2 -> FixParameter(1,87.6488);
    fun2 -> FixParameter(2,5.12017);
    //Integration in three parts
    A = fun1 -> Integral(50,l);
    B = fun1 -> Integral(u,130);
    C = fun1 -> Integral(l,u);
    area1 = C;
    area2 = A+B;
    factor = area1/area2;
    cout << endl;
    cout <<" A: " << A << endl;
    cout <<" B: " << B << endl;
    cout <<" C: " << C << endl;
    cout <<" area1: " << area1 << endl;
    cout <<" area2: " << area2 << endl;
    cout <<" factor: " << factor << endl;
    //To normalize the "pt_BG" in the scale of the background in "pt_signal"
    //pt_BG -> Scale(factor);
    TH1F *pt_BG_rescale = new TH1F("pT_BG Rescale","pT_BG Rescale",50,0,50);  //Change the cutting number of binW from 500 into 100
    pt_BG_rescale -> Add(pt_BG,factor);
    //The subtraction of two graphs
    TH1F *pt_pure = (TH1F*)pt_signal->Clone("pt_pure");
    pt_pure -> Add(pt_BG_rescale,-1.);
    pt_pure -> SetTitle("P_{T} pure");
    //Get some important statistical variables in pt_pure
    cout <<" Entries_pure: " <<  pt_pure -> GetEntries()<< endl;
    cout <<" Mean_pure: " <<  pt_pure -> GetMean()<< endl;
    cout <<" StdDev_pure: " <<  pt_pure -> GetStdDev()<< endl;
    cout <<" BinWidth_pure: " <<  pt_pure -> GetBinWidth(0)<< endl;
    cout <<" Entries_signal: " <<  pt_signal -> GetEntries()<< endl;
    cout <<" Mean_signal: " <<  pt_signal -> GetMean()<< endl;
    cout <<" StdDev_signal: " <<  pt_signal -> GetStdDev()<< endl;
    cout <<" Entries_BG_rescale: " <<  pt_BG_rescale -> GetEntries()<< endl;   //Same as pt_BG
    cout <<" Mean_BG_rescale: " <<  pt_BG_rescale -> GetMean()<< endl;
    cout <<" StdDev_BG_rescale: " <<  pt_BG_rescale -> GetStdDev()<< endl;
    //Some draw options
    pt_BG_rescale -> SetLineColor(kBlack);
    pt_BG_rescale -> SetMarkerColor(kBlack);
    pt_BG_rescale -> SetMarkerStyle(23);
    pt_pure -> SetLineColor(kRed+2);
    pt_pure -> SetMarkerColor(kRed+2);
    pt_pure -> SetMarkerStyle(20);
    //Draw the histogram
    TH1F *pt_distributions = (TH1F*)pt_signal->Clone("pt_distributions");
    pt_distributions -> SetTitle("p_{T} Distributions");
    pt_signal -> SetTitle("p_{T} Signal");
    pt_BG -> SetTitle("p_{T} BG");
    pt_BG_rescale -> SetTitle("p_{T} BG rescale");
    TCanvas *f = new TCanvas("Results_final","Results_final",600,600);
    f -> SetGrid();
    pt_distributions -> Draw();
    pt_BG_rescale -> Draw("SAME");
    pt_pure -> Draw("SAME");   //Whether to fitting or not?
    //Write the graphs into "Results.root"
    TFile *output = new TFile("Results_final.root","recreate");
    pt_signal -> Write();
    pt_BG -> Write();
    pt_BG_rescale -> Write();
    pt_pure -> Write();
	//Add the axes
	TAxis *Xaxis = pt_distributions -> GetXaxis();
	TAxis *Yaxis = pt_distributions -> GetYaxis();
    Xaxis -> SetTitle("p_{T}[GeV/c]");Xaxis -> SetTitleOffset(1.5);Xaxis -> CenterTitle(true);
    Yaxis -> SetTitle("Number of events");Yaxis -> SetTitleOffset(1.5);Yaxis -> CenterTitle(true);
    pt_distributions -> GetYaxis() -> SetRangeUser(-1000,50000);   //Set the Y range
    //pt_distributions -> SetAxisRange(-1000,50000,"Y");
    //Yaxis -> SetRange(-1000,50000);
	//Draw the legend
	TLegend *legend2=new TLegend(0.55,0.65,0.85,0.85);
	legend2->SetTextFont(42);
	legend2->SetTextSize(0.04);
	legend2->AddEntry(pt_signal,"p_{T} Signal+Background","lpe");  //"l" means line;"C" means center;"lpe" or "lep" is used in data with error bar
	legend2->AddEntry(pt_BG_rescale,"p_{T} Background Rescale","lpe");
    legend2->AddEntry(pt_pure,"p_{T} Pure signal","lpe");
	legend2->Draw("same");
    TLegend *legend3=new TLegend(0.55,0.4,0.85,0.6);
    legend3->SetTextFont(42);
    legend3->SetTextSize(0.04);
    legend3->SetHeader("p_{T} Pure","C");
    legend3->AddEntry((TObject*)0,Form("#mu = %.3f #pm %.3f",pt_pure -> GetMean(),pt_pure -> GetMeanError()),"");
    legend3->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f",pt_pure -> GetStdDev(),pt_pure -> GetStdDevError()),"");
    legend3->AddEntry((TObject*)0,Form("BinW = %.1f",pt_pure -> GetBinWidth(0)),"");
    legend3->Draw("same");
}
/*void writehistlist()
{
TList *l = new TList();
TFile *file= new TFile("Results_final.root");           //Results.root is the name of the file.
TH1F *pt_signal = (TH1F*)file -> Get("pt_signal");   //redirect to the graph in "Results.root"
TH1F *pt_BG = (TH1F*)file -> Get("pt_BG");
TH1F *pt_pure = (TH1F*)file -> Get("pt_pure");
l->Add(pt_signal);
l->Add(pt_BG);
l->Add(pt_pure);
TFile *output = new TFile("Results_final.root","RECREATE");
l->Write("Results_final", TObject::kSingleKey);
}*/
