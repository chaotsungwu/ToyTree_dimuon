#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TArray.h>
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
void toytree_final_another_fit()   //we always name the main function the name same as file
{  
    gStyle->SetOptStat(0);   //strip off the original numerical legend.
    TH1F *h_mass[50];
    TCanvas *c[50];
    TLegend *legend1[50];
    TAxis *Xaxis1[50];
    TAxis *Yaxis1[50];
    for (int j=0;j<50;j++){
	    TFile *file = new TFile("Results_final_another.root","read");   //ToyTree.root is the name of the file.
	    h_mass[j] = (TH1F*)file -> Get(Form("pT from %d~%d",j,j+1));   //redirect to the graph of Mass Reconstruction in "Results_final_another.root"
	    h_mass[j] -> Sumw2();   //To tranfer the data of histogram into the data points
    }
    //Create an empty array to fill data in
    double data_x[50];
    double data_x_error[50];
    double data_y[50];
    double data_y_error[50];
    double par[5];
    double ndf = 0.0;
    double chi_squr = 0.0;
    int bin80;   //Find the bin number of x=80.
    double content;
    double error;
    double chi_ndf = 0;
    double error_2 = 0;
    double error_3 = 0;
    double error_4 = 0;
    //因為50張圖都要用同樣的函數去fit，所以只要定義一個function就好，到時候如果真的有幾張fit不好，大不了再定義函數去個別fit，總而言之，不用一開始就定義50個函數。
    //fitting function
    TF1 *fun = new TF1("fun",fitfunction,50,140,5);   //("name",function,function down bound,function up bound,number of parameters)
    //Set the limit of parameter
    fun->SetParLimits(0,10,500000);   //(par0,db,un)
    fun->SetParLimits(1,-1,-0.001);   //(par1,db,un)
    fun->SetParLimits(2,100,500000);   //(par2,db,un)
    fun->SetParLimits(3,80,95);   //(par3,db,un)
    fun->SetParLimits(4,1,6);   //(par4,db,un)
    TF1 *fun1 = new TF1("fun1",GAU,50,130,3);
    TF1 *fun2 = new TF1("fun2",EXP,50,130,2);
    TFile *output = new TFile("Results_final_another_fit.root","recreate");
    for(int j=7;j<36;j++)
    {
	    c[j] = new TCanvas(Form("pT from %0.1d~%0.1d",j,j+1),"",600,600);
        gStyle->SetLegendBorderSize(0);   //Strip off the legend box
        //Using error/content value of one of bins to determine rebin or not.
	    bin80   = h_mass[j] -> GetXaxis() -> FindBin(80);   //Find the bin number of x=80.
	    content = h_mass[j] -> GetBinContent(bin80);
	    error   = h_mass[j] -> GetBinError(bin80);
	    if((error/content)>0.3) h_mass[j] -> Rebin(2);
	    binW = h_mass[j] -> GetBinWidth(0);
        //Fitting each graphs
	    h_mass[j] -> Fit("fun","","",50,130);
	    h_mass[j] -> Draw();
	    //get "Mean" and "chi_squr/ndf"
	    fun -> GetParameters(par);
	    chi_squr = fun -> GetChisquare();
	    ndf = fun -> GetNDF();
	    chi_ndf = chi_squr/ndf;
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
	    //some setting for plot
	    h_mass[j] -> SetTitle(Form("p_{T} : %d #leq p_{T} < %d(Gev/c);mass[Gev/c^{2}];Number of events",j,j+1));
	    TLegend *lg = new TLegend(0.5,0.4,0.75,0.85);   //(X left,Y down,X right,Y up)
        //lg->SetHeader("Some informations","C");
        lg->SetTextSize(0.03);
        //lg->SetTextFont(22);
	    lg->AddEntry(h_mass[j],"data","lpe");
	    lg->AddEntry(fun,"Fitting line","l");
        lg->AddEntry(fun1,"Pure signal","f");
        lg->AddEntry(fun2,"Background","l");
        lg->AddEntry((TObject*)0,Form("N = %.3f #pm %.3f",par[2],error_2),"");
        lg->AddEntry((TObject*)0,Form("M_{#mu#mu} = %.3f #pm %.3f",par[3],error_3),"");
        lg->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f",par[4],error_4),"");
	    lg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.3f",chi_ndf),"");
        lg->AddEntry((TObject*)0,Form("BinW = %.1f",binW),"");
	    lg->Draw("same");
	    //Write the graphs into "Results_final_another_fit.root"
	    c[j] -> Write();
        //c[j] -> Print(Form("c%d.png",j),"png");
	    //上面打lg的原因只是方便你在迴圈裡面把所有圖存在電腦裡面:用 c1 -> SaveAs("地址/名稱");
	    //print on the screen
	    cout << endl;
	    cout <<" Number of events: " << par[2] << endl;//是面積嗎？要想清楚唷～
	    cout <<" Mean: " << par[3] << endl;
	    cout <<" chi_squr/ndf: " << (double) chi_squr/ndf << endl;
        //cout << " binW: " << binW << endl;
	    //Fill the # of events into the array
	    data_x[j] = ((j)+(j+1.0))/2;
	    data_y[j] = par[2];
	    data_y_error[j] = fun -> GetParError(2);
	    data_x_error[j] = 0.5;   //你取pT區間的中間點當作的該區間（寬度為1）代表，那一個點要表示那個區域就必須加上的"誤差"0.5（在這裡其實也不是誤差啦，只是找不到精確的字眼去表示）
        //cout <<"x: "<< data_x[j]<< "     " << "y: " << data_y[j]<<endl;
    }
    TGraphErrors *gpt = new TGraphErrors (36,data_x,data_y,data_x_error,data_y_error);   //Draw the data calculate from above in the graph
    gpt->SetMarkerStyle(8);
    gpt->SetMarkerColor(1);
    gpt->SetTitle("p_{T} distribution;p_{T}[GeV/c];Number of events");
    gpt -> GetXaxis() -> SetRangeUser(5,40);   //Set the X range
    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    //c2->cd();
    gpt->Draw("AP");
    //Write the graphs into "Results_final_another_fit.root"
    gpt -> Write();
    c2 -> Write();
    TCanvas *c3 = (TCanvas*)c2-> Clone("c3");
    c3 -> SetGrid();
    c3 -> Write();
    //C++不像python一樣可以這樣做，要印出array的內容，要寫迴圈，不然就是印出array的開頭位址而已
    //cout << " data_x: " << data_x << endl;
    //cout << " data_y: " << data_y << endl;
    //It should be like:
    for(int k=0;k<50;k++){   //k<sizeof(data_x)
        cout <<"x: "<< data_x[k]<< "     " << "y: " << data_y[k]<<endl;
    }
    
    //cout<< pt_pure1 -> FindBin(7)<<endl;
    //cout<< pt_pure1 -> FindBin(35)<<endl;
    //Plot the ratio of two ways to draw pT distribution
    TFile *file1 = new TFile("Results_final.root");
    TH1F *pt_pure1 = (TH1F*)file1 -> Get("pt_pure");
    double data_x1[29];
    double data_x_error1[29];
    double data_y_error1[29];
    double ratio[29];
    for(int n=0;n<29;n++){
        ratio[n] = (pt_pure1 -> GetBinContent(n+8))/(data_y[n+7]);
        data_x1[n] = ((7.0+n)+(7.0+n+1.0))/2;
        data_x_error1[n] = 0.5;
        data_y_error1[n] = ratio[n]*sqrt(pow(pt_pure1 -> GetBinError(n+8)/pt_pure1 -> GetBinContent(n+8),2)+pow(data_y_error[n+7]/data_y[n+7],2));
    }
    TGraphErrors *gpt1 = new TGraphErrors (29,data_x1,ratio,data_x_error1,data_y_error1);
    gpt1->SetMarkerStyle(8);
    gpt1->SetMarkerColor(1);
    gpt1->SetTitle("The Ratio of Two Ways to Draw p_{T} Distribution");
    TCanvas *c4 = new TCanvas("c4","c4",600,600);
    gpt1->Draw("AP");
    gPad -> SetLeftMargin(0.13);
    c4 -> SetGrid();
    //Add the axes
    TAxis *Xaxis = gpt1 -> GetXaxis();
    TAxis *Yaxis = gpt1 -> GetYaxis();
    Xaxis -> SetTitle("Mass[GeV/c^{2}]") ; Xaxis -> SetTitleOffset(1.2);
    Yaxis -> SetTitle("Ratio") ; Yaxis -> SetTitleOffset(1.8);
    c4 -> Write();
}
