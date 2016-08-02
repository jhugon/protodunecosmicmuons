//Macro cosmic.C
Double_t cosmicspectrum(Double_t *x,Double_t *par)
{
  Float_t xx =x[0];
  Double_t f = 0.14*1e4*(1/(1+1.1*xx*cos(par[0])/115) + 0.054/(1+1.1*xx*cos(par[0])/850));
  return f;
}

void cosmic()
{
  float theta=0;
  TCanvas *c1 = new TCanvas("c1","Cosmic canvas",800,800);
      c1-> SetLogy(1);
      c1-> SetLogx(1);


    TF1 *fa1 = new TF1("fa1",cosmicspectrum,10,1000,1);
    fa1->SetParameter(0,0);
    //--- theta = 0 --- 
//      TF1 *fa1 = new TF1("fa1","0.14*(1/(1+1.1*x*0.7/115 + 0.054/(1+1.1*x*0.7/850)))",10,1000);
//  TF1 *fa1 = new TF1("fa1","0.14*x**(-2.7)*(1/(1+1.1*x*0.7/115 + 0.054/(1+1.1*x*0.7/850)))",1,1000);

//  TF1 *fa2 = new TF1("fa2","0.14*(1/(1+1.1*x*0.7/115))",10,1000);

    TF1 *fa2 = new TF1("fa2",cosmicspectrum,10,1000,1);
    fa2->SetParameter(0,1.31);
    //--- theta = 1.31 rad = 75 deg --- 
  
    fa1->SetLineColor(2);
    fa1->SetLineStyle(7);

    fa2->SetLineColor(4);
    fa2->SetLineStyle(7);

    fa1->SetTitle("cosmic muon flux at surface");
    fa1->GetXaxis()->SetTitle("E_{#mu} [GeV]");  
    fa1->GetYaxis()->SetTitle("E_{#mu}^{2.7} dN/dE [m^{-2} s^{-1} sr^{-1} (GeV)^{1.7}]");

    fa1->Draw();

    //--- valid only for E_mu > 100/cos(theta) GeV ---
    TF1 *fa1 = new TF1("fa1",cosmicspectrum,100,1000,1);
    fa1->SetLineStyle(1);
    fa1->Draw("SAME");


    fa2->Draw("SAME");

    //--- valid only for E_mu > 100/cos(theta) GeV ---
    TF1 *fa2 = new TF1("fa2",cosmicspectrum,386,1000,1);
    fa2->SetParameter(0,1.31);
    fa2->SetLineColor(4);
    fa2->SetLineStyle(1);
    fa2->Draw("SAME");

   
  
    TLegend *legend=new TLegend (0.7,0.75,0.85,0.85);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->AddEntry(fa1,"#theta = 0^{o}","l");
    legend->AddEntry(fa2,"#theta = 75^{o}","l");
    legend->Draw();

}
