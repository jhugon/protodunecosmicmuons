#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include<fstream>
#include <stdlib.h>
#include <time.h>
#include <iomanip>

void SamplingProgram_10GeV()
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                    PHASE 1: THETA VALUES
 
//The following variables should be initialized as needed by the user
 const int numberofevents = 10000;

 const int maxenergy = 10; //Events only generated up to this energy (GeV)
 //Lower energy range means more accurate energy sampling within that range

 srand(22222); //set seed for random number generation      

 const Double_t w = 210;  //Number of sections to divide ranges into.
                         //More slices means higher resolution, longer runtime.
 //w is a master variable, setting the robustness of the entire program's math


 //the following variables should be left alone
 //=====================================================
 int i = 0; //rand holders
 Double_t p = 0;
 int j = 0;
 int k = 0;
 int kmax = 0;
 int m = 0; //counters
 int n = 0;


 Double_t phi = 0;
 Double_t zenith = 0;//angle samples for text output
 Double_t px = 0;
 Double_t py = 0;
 Double_t pz = 0;//momentum samples for text output
 Double_t kineticenergy = 0;//energy sample for text output
 Double_t x = 0;
 Double_t z = 0;//position samples for text output
 ofstream ophelia; //used to output to the text file

const Double_t y = 141.55;
const Double_t muonmass = 0.105658;
const Double_t pi = TMath::Pi();
const Double_t upper = pi*7/18; //maximum allowed theta value (70 degrees)

 Double_t probtotal = 0; //Used to find normconst.
 Double_t normconst = 0; //probtotal * normconst = RAND_MAX = (2^31 - 1)
 //Normconst guarantees the full range of RNG is used in the double if below.

 Double_t probability[w] = {0}; //Probability of finding any event in each
                                //range of angles.

 int theta[w] = {0}; //Number of events in each theta range stored here.

 Double_t probbase = 0; //Running total of ranges checked for each event.

 TF1 theta_equation("CosSquared","cos(x)*cos(x)",0,upper);

 //setup for the integrator
ROOT::Math::GaussIntegrator ig;
ig.SetRelTolerance(0.001);
ROOT::Math::WrappedTF1 wtheta(theta_equation);
ig.SetFunction(wtheta);
 //======================================================


 cout << "Initiating theta sampling portion of program..." << endl;

//Find true values for probability[], probtotal, and normconst

 while(k < w)
   {
     probability[k] = ig.Integral(k*upper/w,(k+1)*upper/w);
     probtotal = probtotal + probability[k];
     k++;
   }

 normconst = 2147483647/probtotal;

 //theta value sampling
 //==============================
 //==============================
 while(j < numberofevents)
   {
     i = rand();
     k = 0;
     probbase = 0;

     while(k < w)
       {
	 
	 //double if plus an if
	 if(i >= probbase)
	   {
	     probbase = probbase + probability[k]*normconst;
	     if(i <= probbase)
	       {
		 theta[k]++;
		 k = w+1;//safe way of breaking out of the loop
	       }
	     else
	       k++;
	   }
	 else
	   k++;

	 if(k==w)//special case: fits into highest theta[] category
	   {
	     k = k - 1;
	     theta[k]++;
	     k++;
	   }

       }//while(k < w)


     j++;
   }//while(j < numberofevents)
 //==============================
 //==============================

 cout << "Done." << endl;
 sleep(1);
 

 //Dump of all values of theta[w] (optional - for diagnostics)
 /*
 cout << endl << "Incoming information dump for theta" << endl << endl;
 sleep(2);
 cout << "Angle range (degrees)" << "     Number of recorded events" << endl << endl;

 while(m < w)
   {
     cout << setw(8) << m*(70/w) << " to " << setw(8) << (m+1)*(70/w) << ":     " << theta[m] << endl;
     m++;
   }
 */

 //Graphing of all values of theta[w] - also optional (events=1000000,w=210)
 /*
 TCanvas *c1 = new TCanvas("c1","Theta Canvas",800,800);
 TH1D *thetahisto = new TH1D("thetahisto","Theta values vs Cos Squared",w,0,w);
 k = 0;
 while(k < w)
   {
     thetahisto -> Fill(k, theta[k]);
     k++;
   }
 thetahisto -> Draw();
 thetahisto -> GetXaxis() -> SetTitle("cos(x*#pi/540)");
 thetahisto -> GetXaxis() -> CenterTitle();
 thetahisto -> GetYaxis() -> SetTitle("Number of events");
 thetahisto -> GetYaxis() -> CenterTitle();
 TF1 *cos2 = new TF1("cos2","cos(7*pi*x/(210*18))*cos(7*pi*x/(18*210))*7550", 0,w);
 cos2 -> Draw("same");
 */

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //                    PHASE 2: PHI VALUES
//Nah, I'm just messin with ya. It's a flat distribution, but in case I ever
//need it not to be, this placeholder is here.

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //                    PHASE 3: ENERGY VALUES


 //This part of the program became mostly useless, as the bulk of its operations were migrated to phase 4. I keep the code here for posterity and for diagnostic reasons. And because I'm too lazy to migrate the rest of it.




 //Reset/initialize variables as needed.
 i = 0;
 j = 0;
 k = 0;
 probbase = 0;
 probtotal = 0;
 normconst = 0;
 Double_t probrange[w] = {0};
 int energy[w] = {0}; //Guess what goes in here?
 //Ready the integration

TF1 flux("Energy Distribution", ".140*pow(x*(1+3.64/(x*pow(sqrt((cos(([0]*pi)/180)*cos(([0]*pi)/180)+0.102573*0.102573+-0.068287*pow(cos(([0]*pi)/180),0.958633)+0.0407253*pow(cos(([0]*pi)/180),0.817285))/(1+0.102573*0.102573+-0.068287+0.0407253)),1.29))),-2.7)*((1/(1+1.1*x*sqrt((cos(([0]*pi)/180)*cos(([0]*pi)/180)+0.102573*0.102573+-0.068287*pow(cos(([0]*pi)/180),0.958633)+0.0407253*pow(cos(([0]*pi)/180),0.817285))/(1+0.102573*0.102573+-0.068287+0.0407253))/115))+0.054/(1+1.1*x*sqrt((cos(([0]*pi)/180)*cos(([0]*pi)/180)+0.102573*0.102573+-0.068287*pow(cos(([0]*pi)/180),0.958633)+0.0407253*pow(cos(([0]*pi)/180),0.817285))/(1+0.102573*0.102573+-0.068287+0.0407253))/850))", 0, maxenergy);



 //main course
 cout << "Initiating energy sampling portion of program..." << endl;
 k = 0;

     flux.SetParameter(0,n*(70/w));
     ROOT::Math::WrappedTF1 wflux(flux);
     ig.SetFunction(wflux);

 while(k < w) //find true values for probtotal, normconst and probrange[]
   {
     probrange[k] = k*k/(w*w)*maxenergy;
//probtotal = probtotal + ig.Integral(probrange[k],(k+1)*(k+1)/(w*w)*maxenergy);
     k++;
   }//while(k < w)
 /*
 normconst = 2147483647/probtotal;


 while(j < m)
   {
    k = 0;
    i = rand();
    probbase = 0;

    while(k < w)
       {
		 
	 //double if plus an if
	 if(i >= probbase)
	   {
       probbase = probbase + ig.Integral(probrange[k],probrange[k+1])*normconst;
	     if(i <= probbase)
	       {
		 if(k > kmax)
		   kmax = k;
		 energy[k]++;
		 k = w+1;//safe way of breaking out of the loop
	       }
	     else
	       k++;
	   }
	 else
	   k++;

       }//while(k < w)
	 
     j++;
   }//while(j < m)

    n++;
   }//while(n < w)
 */
 cout << "Done." << endl;
 sleep(1);

//Dump of all values of energy[w] (optional - for diagnostics)
/* 
 cout << endl << "Incoming information dump for energy" << endl << endl;
 sleep(2);
 cout << "Energy range (GeV)   " << "     Number of recorded events" << endl << endl;

 m = 0;
 while(m < kmax)
   {
     cout << setw(8) << probrange[m] << " to " << setw(8) << probrange[m+1] << ":     " << energy[m] << endl;
     m++;
   }
*/

//Plotting of all energy values (also optional)
/*
 TCanvas *c2 = new TCanvas("c2","Energy Canvas",800,800);
 TH1D *energyhisto = new TH1D("energyhisto","Energy Ranges",100,0,10);
 k = 0;
 while(k < w)
   {
     energyhisto -> Fill(probrange[k],energy[k]);
     k++;
   }
 energyhisto -> Draw();
*/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //                    PHASE 4: COMMUNICATE WITH LARSOFT

 cout << "Initiating output to text file..." << endl;
//It's actually doing a lot more than this, as you can see in the code below.
//I found it necessary to move part of the energy sampling code to this section

ophelia.open("hamlet_10.txt");
j = 0;
k = 0;
m = 0;
n = 0;

      while(k < w)
	{
	  j = 0;

	  while(j < theta[k])
	    {
	      m++;

	      p = rand();
              phi = 2*pi*(p/2147483647);
              p = rand();
              zenith =  (k*(70/w)+(70/w)*(p/2147483647))*pi/180;

              px = sin(zenith)*cos(phi);
              pz = sin(zenith)*sin(phi);
              py = cos(zenith); //y is vertical in this detector

	      p = rand();
	      x = 349.666*(p/2147483647) - 78.816;
	      p = rand();
	      z = 231.065*(p/2147483647) - 38.77;

              flux.SetParameter(0,k*(70/w)+(35/w));
              ROOT::Math::WrappedTF1 wflux(flux);
              ig.SetFunction(wflux);
	      n = 0;
	      probtotal = 0;

	      while(n < (w-1))//The code will break if it ever gets to n = w.
		{         //For appropriate values of maxenergy, this will
                          //never ever happen.
	       probtotal = probtotal + ig.Integral(probrange[n],probrange[n+1]);
	       //Take note of the fact that probrange[] covers progressively
	       //larger energy ranges. This is because the integral function
	       //loses a bit of accuracy over larger ranges, so the lower
	       //energy ranges (which contain the vast majority of events, and
	       //are more important to us) have very small ranges to maintain
	       //accuracy, and the higher energies are lumped into larger
	       //ranges because so few of those events will be generated, and
	       //aren't as important contextually for the near detector
	       n++;
		}

	      normconst = 2147483647/probtotal;
	      n = 0;
	      probbase = 0;
	      i = rand();

	      while(n < w)
		{

		  //double if
		  if(i >= probbase)
		    {
       probbase = probbase + ig.Integral(probrange[n],probrange[n+1])*normconst;
		      if(i <= probbase)
			{
			  p = rand();
      kineticenergy = probrange[n]+(probrange[n+1]-probrange[n])*(p/2147483647);
			  n = w+1;//safe way of breaking out of the loop
			}
	               else
	               n++;
		    }
	          else
		  n++;

		}//while(n < w)

kineticenergy = TMath::Power(kineticenergy*kineticenergy+muonmass*muonmass,0.5);
//Event number and number of particles in event (always 1)
 ophelia << m << " 1" << endl;
      //Status code (1 so that particle is tracked) and PDG code of particle (13 for muon, -13 for antimuon)
 //0's to indicate no parent particles or daughter particles
 ophelia << "1 -13 0 0 0 0 ";
   //The components of the 4-momentum (px, py, pz, t) and the particle rest mass
  ophelia  <<  px*kineticenergy << " " << py*kineticenergy << " "  << pz*kineticenergy << " " << kineticenergy << " " << muonmass << " ";
     
//Position and time coordinates of the particle's starting position (x, y, z, t)
  ophelia << x << " " << y << " " << z << " 0.0" << endl;	    
 j++;
	    }//while(j < theta[k])
	  k++;
	}//while(k < w)

ophelia.close();
cout << "Done." << endl;

}//void SamplingSally(void)
