//This macro uses a modified version of Gaisser's formula, found in the article
//"Muon simulation at the Daya Bay site", which can be found at the url
// http://escholarship.org/uc/item/6jm8g76d
//The modified formula includes a term to account for the possibility of
//muons decaying during their flight through the atmosphere, and so is able
//to better describe cosmic ray muon flux at sea level for lower energy muons


Double_t fluxformula(Double_t *x, Double_t *par)
{
  //x = energy of the muons, par = zenith angle
  Double_t xx = x[0];
  Double_t p1 = 0.102573;
  Double_t p2 = -0.068287;
  Double_t p3 = 0.958633;
  Double_t p4 = 0.0407253;
  Double_t p5 = 0.817285;
  Double_t costheta = cos((par[0]*TMath::Pi)/180);
  Double_t workaround = TMath::Power(costheta,p5);
  //the term on the above line was acting up when put into the parstar formula

  Double_t parstar = TMath::Sqrt((costheta*costheta+p1*p1+p2*TMath::Power(costheta,p3)+p4*workaround)/(1+p1*p1+p2+p4));

  Double_t twopointseven = xx*(1+3.64/(xx*TMath::Power(parstar,1.29)));
  //please excuse my naming conventions, but they help me to remember!

  Double_t f = .140*TMath::Power(twopointseven,-2.7)*((1/(1+1.1*xx*parstar/115))+0.054/(1+1.1*xx*parstar/850));
  return f;
}



void MuonFlux(void)
{
  TCanvas *c1 = new TCanvas("c1","Muon Canvas",500,500);

  //  c1 -> SetLogx(1);
  //c1 -> SetLogy(1); //I saw no need to implement these, but they are available

 TF1 *flux = new TF1("flux",fluxformula,0,10,1);
 flux -> SetParameter(0,0);
 
 TF1 *flux2 = new TF1("flux2",fluxformula,0,10,1);
 flux2 -> SetParameter(0,30);
 flux2 -> SetLineColor(kCyan);

 TF1 *flux3 = new TF1("flux3",fluxformula,0,10,1);
 flux3 -> SetParameter(0,60);
 flux3 -> SetLineColor(kSpring+5);

 flux -> SetTitle("Muon Flux at Sea Level - Modified Formula");
 flux -> GetXaxis() -> SetTitle("E_{#mu} [GeV]");
 flux -> GetYaxis() -> SetTitle("E_{#mu} dN/dE [m^{-2} s^{-1} sr^{-1} GeV]");
 flux -> GetXaxis() -> CenterTitle();
 flux -> GetYaxis() -> CenterTitle();
 flux -> GetYaxis() -> SetTitleOffset(2);
 gStyle -> SetPadLeftMargin(.15);

 flux -> Draw();
 flux2 -> Draw("same");
 flux3 -> Draw("same");

 TLegend *legend = new TLegend (0.7,0.75,0.85,0.85);
 legend -> AddEntry(flux, "#theta = 0^{o}","l");
 legend -> AddEntry(flux2, "#theta = 30^{o}","l");
 legend -> AddEntry(flux3, "#theta = 60^{o}","l");
 legend -> Draw();
}
