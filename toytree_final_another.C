#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TArray.h>
double binW;
void toytree_final_another()   //we always name the main function the name same as file
{  
    gStyle->SetOptStat(0);   //strip off the original numerical legend.
    TFile *file1= new TFile("ToyTree.root");           //ToyTree.root is the name of the file.
    TTreeReader reader("ToyTree",file1);               //ToyTree is the name of tree in the file.
    TTreeReaderArray<double> pt1(reader,"mu1_pt");
    TTreeReaderArray<double> pt2(reader,"mu2_pt");
    TTreeReaderArray<double> eta1(reader,"mu1_eta");
    TTreeReaderArray<double> eta2(reader,"mu2_eta");
    TTreeReaderArray<double> phi1(reader,"mu1_phi");
    TTreeReaderArray<double> phi2(reader,"mu2_phi");
    //分50張圖就夠了
    TH1F *h_mass[50];
    char name[50];
    char title[50];
    for (int j=0;j<50;j++){
        sprintf(name,"pT from %d~%d",j,j+1);
        sprintf(title,"pT from %d~%d",j,j+1);
	//bin width可以大一點
        h_mass[j] = new TH1F(name,title,200,50,150);
        h_mass[j] -> Sumw2();   //To tranfer the data of histogram into the data points
    }
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
    double j=0;
    double ndf = 0.0;
    double chi_squr = 0.0;
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
            for(int j=0;j<50;j++)
            {
                if (pt > j && pt < (j+1) ){
                    h_mass[j] -> Fill(Mass);
                }
            }
        }
    }
    //Draw the histogram
    TFile *output= new TFile("Results_final_another.root","recreate");
    for(int j=0;j<50;j++){
           h_mass[j] -> Write();
    }
}
