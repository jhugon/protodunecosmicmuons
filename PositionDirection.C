//put these somewhere near the top of the code - probably next to all the other
//variables in the 'leave these alone'  section
//======================================================
Double_t phi = 0;
Double_t zenith = 0; //incoming angles, used to find the initial direction

Double_t x = 0;
Double_t y = 0;
Double_t z = 0; //Starting position of the particle

Double_t px = 0;
double_t py = 0;
Double_t pz = 0; //Initial momentum; dependent on initial direction and energy

Double_t kineticenergy = 0; //LArSoft needs a specific energy, not a range
const Double_t electronmass = 0.511; //Rest energy of an electron, in MeV
const Double_t pi = TMath::Pi();
Double_t totalenergy = 0;
//=====================================================


//Put this right below the line that says 'while(j < numberofevents)'
//But not before the { on the line immediately after the while!
//{s are important and should always immediately follow the if() or while() line
//So it should look like this:

//while(j < numberofevents)
//{
//[The code you insert]


//=====================================================
p = rand();
phi = 2*pi*(p/2147483647); //Random phi angle from 0 to 2pi
p = rand();
zenith = pi*(p/2147483647);//Random zenith angle from 0 to pi

px = sin(zenith)*cos(phi);
pz = sin(zenith)*sin(phi);
py = cos(zenith); //y is vertical in LarSoft
//sqrt(px^2 + py^2 + pz^2) should equal 1 at this point

p = rand();
x = 294.666*(p/2147483647) - 51.316; //random x value between xmin and xmax
p = rand();
z = 176.065*(p/2147483647) - 11.27; //random z value between zmin and zmax
p = rand();
y = 228.07*(p/2147483647) - 96.52; //random y value between ymin and ymax
//=====================================================




//The next bit of code needs to go in two places
//One is smack dab in the middle of the double if statement
//The other is within the third if statement below those two
//It should look like this
/*
if(i >= probbase)
	   {
	     probbase = probbase + probability[k]*normconst;
	     if(i <= probbase)
	       {
		 events[k]++;
FIRST PLACE TO INSERT THE CODE
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
SECOND PLACE TO INSERT THE CODE
	     k++;
	   }
*/


//=====================================================
p = rand();
kineticenergy = k*(w/45) + (w/45)*(p/2147483647);
//The above line picks a random energy from the currently selected range
//For example, if the event randomly fell into the energy range between
//6 and 7 MeV (for w=45), the above code would numerically look like this
//kineticenergy = 6*(45/45) +(45/45)*[random number between 0 and 1]
//kineticenergy = 6*1 + 1*[random number between 0 and 1]
//kineticenergy = 6 + [random number between 0 and 1]
//So kineticenergy will be given a specific energy between 6 and 7 MeV
//You don't have to include all these comments in the code unless you want to
//=====================================================



//And finally, the output section, with a few finishing touches on the variables
//Put this directly above the j++ line, like so:

/*
	 if(k==w)//special case: fits into highest events[] category
	   {
	     k = k - 1;
	     events[k]++;
	     k++;
	   }

       }//while(k < w)

INSERT THE OUTPUT CODE HERE
     j++;
   }//while(j < numberofevents)
 //==============================
 cout << "Done." << endl;
 sleep(1);
*/


totalenergy = TMath::Power(kineticenergy*kineticenergy+electronmass*electronmass,0.5);
//Event number and number of particles in event (always 1)
 ophelia << j+1 << " 1" << endl;
//Status code (1 so that particle is tracked) and PDG code of particle (11 for electron)
 //0's to indicate no parent particles or daughter particles
 ophelia << "1 11 0 0 0 0 ";
//The components of the 4-momentum (px, py, pz, t) and the particle rest mass
  ophelia  <<  px*kineticenergy << " " << py*kineticenergy << " "  << pz*kineticenergy << " " << kineticenergy << " " << electronmass << " ";
//Position and time coordinates of the particle's starting position (x, y, z, t)
  ophelia << x << " " << y << " " << z << " 0.0" << endl;



//P.S. you can delete the code that's commented out at the bottom of the program that looks similar to the above stuff
//P.P.S. put this code on the second to last line, right above the final }

ophelia.close();
