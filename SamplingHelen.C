#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include<fstream>
#include <stdlib.h>
#include <time.h>
#include <iomanip>

void SamplingHelen()
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                    PHASE 1: THETA VALUES
 
//The following variables should be initialized as needed by the user
 const int numberofevents = 1000;

 srand(19850); //set seed for random number generation      

 const Double_t w = 46;  //Number of sections to divide ranges into.
                         //More slices means higher resolution, longer runtime.
 //w is a master variable, setting the robustness of the entire program's math


 //the following variables should be left alone
 //=====================================================
 int i = 0; //rand holders
 Double_t p = 0;
 int j = 0;
 int k = 0;

 ofstream ophelia; //used to output to the text file

 Double_t probtotal = 0; //Used to find normconst.
 Double_t normconst = 0; //probtotal * normconst = RAND_MAX = (2^31 - 1)
 //Normconst guarantees the full range of RNG is used in the double if below.

 Double_t probability[w] = {0}; //Probability of finding any event in each
                                //range of energies.

 int events[w] = {0}; //Number of events in each energy range stored here.

 Double_t probbase = 0; //Running total of ranges checked for each event.

 TF1 neutrino_energy_spectrum("Neutrino Energy Spectrum","x*x",0,45); //Replace x*x with your own equation

//setup for the integrator
ROOT::Math::GaussIntegrator ig;
ig.SetRelTolerance(0.001);
ROOT::Math::WrappedTF1 wnes(neutrino_energy_spectrum);
ig.SetFunction(wnes);
 //======================================================


 cout << "Initiating sampling portion of program..." << endl;

//Find true values for probability[], probtotal, and normconst

 while(k < w)
   {
     probability[k] = ig.Integral(k*(45/w),(k+1)*(45/w));
     probtotal = probtotal + probability[k];
     k++;
   }

 normconst = 2147483647/probtotal;

 //energy value sampling
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
		 events[k]++;
		 k = w+1;//safe way of breaking out of the loop
	       }
	     else
	       k++;
	   }
	 else
	   k++;

	 if(k==w)//special case: fits into highest events[] category
	   {
	     k = k - 1;
	     events[k]++;
	     k++;
	   }

       }//while(k < w)


     j++;
   }//while(j < numberofevents)
 //==============================
 cout << "Done." << endl;
 sleep(1);



 //Dump of all values of events[w] (optional - for diagnostics)

 cout << endl << "Incoming information dump" << endl << endl;
 sleep(2);
 
 j = 0;
 cout <<"        Energy Range   Number of events" << endl;
 while(j < w)
   {
    
     cout << setw(8) << j*45/w << " to " << setw(8) << (j+1)*(45/w)  << "   " << events[j] << endl;
     j++;
   }




 /*
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
 */

}//void SamplingSally(void)
