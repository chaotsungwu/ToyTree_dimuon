#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TArray.h>
void toytree_final_combine()   //we always name the main function the name same as file
{  
    gStyle->SetOptStat(0);   //strip off the original numerical legend.
    TFile *file1 = new TFile("Results_final.root");           //Results.root is the name of the file.
    TFile *file2 = new TFile("Results_final_another_fit.root");
    TH1F *pt_pure1 = (TH1F*)file1 -> Get("pt_pure");   //redirect to the graph in "Results_final.root"
    //TGraphErrors *pt_pure2 = (TGraphErrors*)file2 -> Get("Graph");
    TGraphErrors *pt_pure2 = (TGraphErrors*)file2 -> Get("c2");
    //TCanvas *c = new TCanvas("pt_combine","pt_combine",600,600);
    //c -> SetGrid();
    pt_pure2 -> Draw("same");
    //pt_pure2 -> GetXaxis() -> SetRangeUser(0.0001,50);   //Set the X range
    pt_pure1 -> Draw("same");
    pt_pure2 -> SetTitle("P_{T} distribution;P_{T}[GeV/c];Number of events");
    TLegend *legend = new TLegend(0.6,0.65,0.88,0.85);
    legend -> SetTextFont(72);
    legend -> SetTextSize(0.04);
    legend -> SetHeader("Method","C");
    legend -> AddEntry(pt_pure1,"Sideband Subtraction","lpe");
    legend -> AddEntry(pt_pure2,"Fittings","lpe");
    legend -> Draw("same");
    TFile *output = new TFile("Results_final_combine.root","recreate");
    pt_pure2 -> Write();
}
