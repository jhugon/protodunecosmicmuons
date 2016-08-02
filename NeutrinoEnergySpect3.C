//This macro uses a formula cited in Chapter 6 of the LBNE science document, 
//which can be found at "https://sharepoint/fnal/gov/project/lbne/LBNE at Work/
//LBNE Science Program/SitePages/Home.aspx"
//This formula describes the supernova-neutrino spectrum at a given moment
//in time.


Double_t EnergyFormula(Double_t *E, Double_t *par)
{
  //x = energy of the muons, par = zenith angle
  Double_t xx = E[0];
  //Double_t p1 = par[0]; 3
  //Double_t p2 = par[1]; 11
  Double_t f = (TMath::Power((xx/(11)),3))*(TMath::Exp(-(1+3)*(xx/(11))));
  return f;
}



void NeutrinoEnergySpect3(void)
{
  TCanvas *c1 = new TCanvas("c1","Neutrino Canvas",500,500);

  //  c1 -> SetLogx(1);
  //c1 -> SetLogy(1); //I saw no need to implement these, but they are available

  TF1 *energy1 = new TF1("energy1",EnergyFormula,0,45,2);
  //energy1 -> SetParameter(0,3);
  //energy1 -> SetParameter(1,11);
 
 //TF1 *energy2 = new TF1("energy2",EnergyFormula,0,45,2);
 //energy2 -> SetParameter(0,2);
 //energy2 -> SetParameter(1,10);
 //energy2 -> SetLineColor(kCyan);

 //TF1 *energy3 = new TF1("energy3",EnergyFormula,0,45,2);
 //energy3 -> SetParameter(0,2);
 //energy3 -> SetParameter(1,12);
 //energy3 -> SetLineColor(kSpring+5);

 energy1 -> SetTitle("Neutrino Energy Spectrum");
 energy1 -> GetXaxis() -> SetTitle("Energy [MeV]");
 energy1 -> GetYaxis() -> SetTitle("#phi(E) [a.u.] ");
 energy1 -> GetXaxis() -> CenterTitle();
 energy1 -> GetYaxis() -> CenterTitle();
 energy1 -> GetYaxis() -> SetTitleOffset(1.5);

 // energy1 -> SetLineColor(kred);

 energy1 -> Draw();
 //energy2 -> Draw("same");
 //energy3 -> Draw("same");

 // TLegend *legend = new TLegend (0.7,0.75,0.85,0.85);
 // legend -> AddEntry(energy1, "Eavg = 11 MeV, #alpha = 3" ,"l");
 //legend -> AddEntry(energy2, "Eavg = 10 MeV","l");
 //legend -> AddEntry(energy3, "Eavg = 12 MeV","l");
 // legend -> Draw();

 TF1 *energy2 = new TF1("energy2",CSFormula,0,45,2);
 energy2 -> SetParameter(0,2);
 energy2 -> SetParameter(1,8);
 energy2 -> SetLineColor(kSpring);


 energy2 -> Draw("same");
 
 c1->SetLogy(1);

 Float_t rightmax = 5E-42; 
 Float_t rightmin = 8E-45;

 
 // TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),rightmin,rightmax,45,"G");
 // TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),2.2e-2,rightmin,rightmax,450,"G");
 TGaxis *axis = new TGaxis(45,0,45,4e-2,rightmin,rightmax,450,"G");
  axis->SetLineColor(kSpring);
  axis->SetLabelColor(kSpring);
  axis->Draw("same"); 


  //energy2 -> Draw();

}
//This macro uses a formula from arxiv.org/pdf/hep-ph/0307222v1.pdf
//This formula describes the cross-section of the neutrino spectrum (this formula does not contain the F(Ee) term) 


Double_t CSFormula(Double_t *E, Double_t *par)
{
  Double_t xx = E[0]; //Energy of the neutrinos
  Double_t p1 = par[0]; //alpha, the pinching parameter
  Double_t p2 = par[1]; //Average energy of the neutrinos
  Double_t p3 = xx-5.885+0.510998910; //Ee term in the formula, =Eneutrino-Q+me
  Double_t p4 = 0.510998910; //me, the mass of the electron 
  //  Double_t f = (1.72)*(TMath::Power(10,-44))*(p3/(p4*p4))*(TMath::Power(((p3*p3)-(p4*p4)),(1/2)));
  Double_t f = 4.e-2/5e-42*(1.72)*(TMath::Power(10,-44))*(p3/(p4*p4))*(TMath::Power(((p3*p3)-(p4*p4)),(1/2)));
  //Double_t f = sin(xx)*p1*p2*p3*p4; This term was used to test the parameters and locate an error
  return f;
}



void CrossSectionTest3(void)
{
  cout << "CJ was here" << endl;
  // TCanvas *c1 = new TCanvas("c1","Neutrino Canvas",500,500);
  
  //c1 -> SetLogx(1);
  // c1 -> SetLogy(1);  

 TF1 *energy1 = new TF1("energy1",CSFormula,0,45,2);
 energy1 -> SetParameter(0,2);
 energy1 -> SetParameter(1,8);
 energy1 -> SetLineColor(kSpring);
 
 
 //TF1 *energy2 = new TF1("energy2",CSFormula,0,45,2);
 //energy2 -> SetParameter(0,2);
 //energy2 -> SetParameter(1,10);
 //energy2 -> SetLineColor(kCyan);

 
 //TF1 *energy3 = new TF1("energy3",CSFormula,0,45,2);
 // energy3 -> SetParameter(0,2);
 // energy3 -> SetParameter(1,12);
 // energy3 -> SetLineColor(kSpring);

 
 //energy1 -> SetTitle("Neutrino Cross Section");
 //energy1 -> GetXaxis() -> SetTitle("Energy [MeV]");
 //energy1 -> GetYaxis() -> SetTitle("#sigma(E) [a.u.] ");
 //energy1 -> GetXaxis() -> CenterTitle();
 //energy1 -> GetYaxis() -> CenterTitle();
 //energy1 -> GetYaxis() -> SetTitleOffset(1.5);
 //energy1 -> Draw("same");
 //c1->Update();
 //energy2 -> Draw("same");
 //energy3 -> Draw("same");

 //TLegend *legend = new TLegend (0.7,0.75,0.85,0.85);
 //legend -> AddEntry(energy1, "Eavg = 8 MeV" ,"l");
 //legend -> AddEntry(energy2, "Eavg = 10 MeV","l");
 //legend -> AddEntry(energy3, "Eavg = 12 MeV","l");
 //legend -> Draw();

 Float_t rightmax = 1E-42; //1.1*energy1->GetMaximum();
 Float_t rightmin = 1E-45;
 Float_t scale = gPad->GetUymax()/rightmax;
 energy1 -> Scale(scale);
 energy1 -> Draw("same");
 
 TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),rightmin,rightmax,45,"+L");
 axis->SetLineColor(kSpring);
 axis->SetLabelColor(kSpring);
 c1->Update();
 axis->Draw("same");

}
