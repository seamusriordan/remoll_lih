#include "remollGenMott.hh"

#include "CLHEP/Random/RandFlat.h"

#include "Randomize.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "remollBeamTarget.hh"
#include "remollMultScatt.hh"
#include "remolltypes.hh"

#include <math.h>

#define Euler 0.5772157
#define NINTERVAL 3

remollGenMott::remollGenMott(){
    fTh_min =     3.*deg;
    fTh_max =     6.*deg;

    fE_min = 1.0*GeV; // Absolute minimum of electron energy
                            // to generate

    fApplyMultScatt = true;
    fBeamTarg = remollBeamTarget::GetBeamTarget();
}

remollGenMott::~remollGenMott(){
}

void remollGenMott::SamplePhysics(remollVertex *vert, remollEvent *evt){
    // Generate ep event
    
    //  Crazy weighting for brem because ep cross section blows up at low Q2

    // Get initial beam energy instead of using other sampling
    double beamE = fBeamTarg->fBeamE;
    double Ekin  = beamE - electron_mass_c2;

    std::vector <G4VPhysicalVolume *> targVols = fBeamTarg->GetTargVols();

    bool bypass_target = false;


    double bremcut = fBeamTarg->fEcut;

    // Approximation for Q2, just needs to be order of magnitude
    double effQ2 = 2.0*beamE*beamE*(1.0-cos(0.5*deg));

    // About ~1.5%
    double int_bt = 0.75*(alpha/pi)*( log( effQ2/(electron_mass_c2*electron_mass_c2) ) - 1.0 );

    double bt;
    double this_radlen = 1e9;

    double ext_len = fBeamTarg->fTravLen/this_radlen;

    bt = (4.0/3.0)*(ext_len + int_bt);

    double prob, prob_sample, sample, eloss, value;
    value = 1.0;
    prob = 1.- pow(bremcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(bremcut/Ekin,bt+1.))
	+ 0.75*bt/(2.+bt)*(1.- pow(bremcut/Ekin,bt+2.));
    prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */
    prob_sample = G4UniformRand();        /* Random sampling */

    double Evlo[NINTERVAL] = {
	bremcut,
	(beamE-bremcut)*2.0*GeV/(11.0*GeV-bremcut),
	(beamE-bremcut)*9.0*GeV/(11.0*GeV-bremcut),
    };

    double Evhi[NINTERVAL] = {
	(beamE-bremcut)*2.0*GeV/(11.0*GeV-bremcut),
	(beamE-bremcut)*9.0*GeV/(11.0*GeV-bremcut),
	(beamE-bremcut)*(11.0*GeV-fE_min)/(11.0*GeV-bremcut),
    };

    assert( Evhi[NINTERVAL-1]-Evlo[NINTERVAL-1] > 0.0 );

    double Eprob[NINTERVAL]  = { 0.40, 0.20, 0.40 };

    double Enorm[NINTERVAL];
    // Interval normalization
    for( int idx = 0; idx < NINTERVAL; idx++ ){
	Enorm[idx]  = ((Evhi[idx]-Evlo[idx])/(Evhi[NINTERVAL-1]-Evlo[0]))
	    /Eprob[idx];
    }

    int    Evidx;
    double evsum = 0.0;
    double vweight = 0.0;
    eloss = 0.0;
	     
    // Averages over the intervals
    double vavg[NINTERVAL] = {
	log(Evhi[0]/Evlo[0])/(Evhi[0]-Evlo[0]),
	log((beamE-Evlo[1])/(beamE-Evhi[1]))/(Evhi[1]-Evlo[1]),
	(1.0/(beamE-Evhi[2])-1.0/(beamE-Evlo[2]))/(Evhi[2]-Evlo[2])
    };

    if (prob_sample <= prob) {//Bremsstrahlung has taken place!
	//  We break this into 4 seperate energy loss intervals
	//  with total integrals roughly the size of
	//  what the ep product looks like with 11 GeV beam
	//   cut  -  2000 MeV, 1/x, 40%
	//   2000 -  9000 MeV, 1/(E-x), 20%
	//   9000 - 10990 MeV, 1/(E-x)^2, 40%

	sample = G4UniformRand();

	// Identify our region
	// based on the probability distribution
	Evidx = 0;
	evsum  = Eprob[Evidx];
	while( evsum < sample ){
	    Evidx++;
	    evsum += Eprob[Evidx];
	}

	sample = G4UniformRand();

	if( Evidx == 0 ){
	    eloss = Evlo[Evidx]*pow(Evhi[Evidx]/Evlo[Evidx],sample);
	    vweight = eloss;
	}

	if( Evidx == 1 ){
	    eloss = beamE - (beamE-Evhi[Evidx])*
		pow((beamE-Evlo[Evidx])/(beamE-Evhi[Evidx]),sample);
	    vweight = (beamE-eloss);
	}

	if( Evidx == 2 ){
	    eloss = beamE - pow( (1.0/(beamE-Evhi[Evidx]) - 1.0/(beamE-Evlo[Evidx]))*sample
		    + 1.0/(beamE-Evlo[Evidx]), -1.0 );
	    vweight = (beamE-eloss)*(beamE-eloss);
	}

	if( !(eloss > 0.0 ) ){
	    printf("idx = %d\n", Evidx );
	}

	assert( eloss > 0.0 );

	assert( !std::isnan(eloss) && !std::isinf(eloss) );

	vweight *= vavg[Evidx];
	//  mult by beamE-bremcut for proper normalization
	value = RadProfile( eloss, bt)*
	    ((Evhi[NINTERVAL-1]-Evlo[0])/EnergNumInt(bt, Evlo[0], Evhi[NINTERVAL-1])) // average of RadProfile
	    *vweight // sampling weighting (flat or ) / average value for normalization
	    *Enorm[Evidx]; //  Weight given the region

	beamE -= eloss;
    }

    if( beamE < electron_mass_c2 ){ 
	evt->SetEffCrossSection(0.0);
	evt->SetAsymmetry(0.0);
	return; 
    }

    // Set event information to our new sampling
    evt->fBeamE = beamE;
    evt->fBeamMomentum = evt->fBeamMomentum.unit()*sqrt(beamE*beamE - electron_mass_c2*electron_mass_c2);;

    ////////////////////////////////////////////////////////////////////////////////////////////


    // sample with 1.0/(1-cos)^2

    double cthmin = cos(fTh_min);
    double cthmax = cos(fTh_max);

    double icth_b = 1.0/(1.0-cthmax);
    double icth_a = 1.0/(1.0-cthmin);

    double sampv = 1.0/CLHEP::RandFlat::shoot(icth_b, icth_a);

    assert( -1.0 < sampv && sampv < 1.0 );

    double th = acos(1.0-sampv);

    // Value to reweight cross section by to account for non-uniform
    // sampling
    double samp_fact = sampv*sampv*(icth_a-icth_b)/(cthmin-cthmax);

    double ph = CLHEP::RandFlat::shoot(0.0, 2.0*pi);

    double ef    = proton_mass_c2*beamE/(proton_mass_c2 + beamE*(1.0-cos(th)));;

    double q2  = 2.0*beamE*ef*(1.0-cos(th));
    double tau = q2/(4.0*proton_mass_c2*proton_mass_c2);

    double gd = pow( 1.0 + q2/(0.71*GeV*GeV), -2.0 );
    double gep = gd;
    double gmp = 2.79*gd;

    double gen =  1.91*gd*tau/(1.0+5.6*tau); // galster
    double gmn = -1.91*gd;

    double sigma_mott = hbarc*hbarc*pow(alpha*cos(th/2.0), 2.0)/pow(2.0*beamE*sin(th/2.0)*sin(th/2.0), 2.0);
    double ffpart1 = (gep*gep + tau*gmp*gmp)/(1.0+tau);
    double ffpart2 = 2.0*tau*gmp*gmp*tan(th/2.0)*tan(th/2.0);

    double sigma_proton = sigma_mott*(ef/beamE)*(ffpart1 + ffpart2);

    // Get total number of atoms

    // Randomly sample
    int thiselidx = CLHEP::RandFlat::shootInt(vert->GetMaterial()->GetNumberOfElements());
    const G4Element *thisel = vert->GetMaterial()->GetElement(thiselidx);

    // FIXME:  Put in A-dependent cross section
    // with thisel->GetA()

    // Annu Rev Nucl Sci 1957.7:231-316, eq 129
    double a = pow(thisel->GetA(), 1.0/3.0)*1.2*fermi;
    double alpha_N = (thisel->GetZ() - 2.0)/3.0;
    double a0 = a/(sqrt((3.0*(2.0+5.0*alpha_N))/(2.0*(2.0+3.0*alpha_N))));
    double q = sqrt(q2)/hbarc;
    double F_generic = (1.0 - ( ( alpha_N*q*q*a0*a0)/(2.0*(2.0+3.0*alpha))))*exp(-q*q*a0*a0/4.0); 

    double sigma = 0.0;

    if( thisel->GetZ() == 1 && thisel->GetA() < 1.5 ){
        sigma = sigma_proton;
    } else {
        sigma = sigma_mott*(ef/beamE)*F_generic*F_generic;
    }

    double V = 2.0*pi*(cthmin - cthmax)*samp_fact;

    // Suppress too low angles from being generated
    // If we're in the multiple-scattering regime
    // the cross sections are senseless.  We'll define this 
    // as anything less than three sigma of the characteristic
    // width
    
    if( th < 3.0*fBeamTarg->fMS->GetPDGTh() ){
	sigma = 0.0;
    }

    //  Multiply by Z because we have Z protons 
    //  value for uneven weighting

    double thisZ = thisel->GetZ();

    // Weight by number of atoms as effective luminosity is in molecules
    evt->SetEffCrossSection(sigma*V*thisZ*thisZ*value*vert->GetMaterial()->GetAtomsVector()[thiselidx]);

    G4double APV_base = -GF*q2/(4.0*sqrt(2.0)*pi*alpha);

    G4double eps = pow(1.0 + 2.0*(1.0+tau)*tan(th/2.0)*tan(th/2.0), -1.0);

    G4double apvffnum = eps*gep*gen + tau*gmp*gmn;
    G4double apvffden = eps*gep*gep  + tau*gmp*gmp;

    G4double APV = APV_base*(QWp - apvffnum/apvffden);

    evt->SetAsymmetry(APV);

    evt->SetQ2( q2 );
    evt->SetW2( proton_mass_c2*proton_mass_c2 );

    // REradiate////////////////////////////////////////////////////////////////////////////
    // We're going to use the new kinematics for this guy

    int_bt = (alpha/pi)*( log( q2/(electron_mass_c2*electron_mass_c2) ) - 1.0 );
    Ekin = ef - electron_mass_c2;;
    double env, ref;

    prob = 1.- pow(bremcut/Ekin, int_bt) - int_bt/(int_bt+1.)*(1.- pow(bremcut/Ekin,int_bt+1.))
	+ 0.75*int_bt/(2.+int_bt)*(1.- pow(bremcut/Ekin,int_bt+2.));
    prob = prob/(1.- int_bt*Euler + int_bt*int_bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */
    prob_sample = G4UniformRand();        /* Random sampling */

    if (prob_sample <= prob) {//Bremsstrahlung has taken place!
	do {
	    sample = G4UniformRand();
	    eloss = fBeamTarg->fEcut*pow(Ekin/fBeamTarg->fEcut,sample);
	    env = 1./eloss;
	    value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);

	    sample = G4UniformRand();
	    ref = value/env;
	} while (sample > ref);

	ef = Ekin-eloss+electron_mass_c2;
	assert( ef > electron_mass_c2 );
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0), 
	    G4ThreeVector(ef*cos(ph)*sin(th), ef*sin(ph)*sin(th), ef*cos(th) ), 
	    "e-" );

    return;

}

G4double remollGenMott::RadProfile(G4double eloss, G4double btt){
     double Ekin = fBeamTarg->fBeamE - electron_mass_c2;
     double retval = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,btt);

     if( std::isnan(retval) || std::isinf(retval) ){
	 G4cerr << __FILE__ << " line " << __LINE__ << ": ERROR" << G4endl;
	 G4cerr << "Ekin " << Ekin/GeV << " GeV   btt = " << btt << " retval = " << retval << G4endl;
	 fprintf(stderr, "eloss = %e GeV\n", eloss/GeV);
     }

     assert( !std::isnan(retval) && !std::isinf(retval) );

     return retval;
}

G4double remollGenMott::EnergNumInt(G4double btt, G4double a0, G4double b0){
    const int nbin = 1000;
    double sum = 0.0;
    double bremcut = fBeamTarg->fEcut;

    int j;
    double boolc[5] = {7.0, 32.0, 12.0, 32.0, 7.0};

    double a, b, thissum;

    for(int i =0;i<nbin;i++) {
	// Integrate over sample spacings that are logarithmic
	a = bremcut*pow(b0/a0, (((double) i)/nbin));
	b = bremcut*pow(b0/a0, (((double) i+1.0)/nbin));

	// Boole's rule
	thissum = 0.0;
	for( j = 0; j < 5; j++ ){
	    thissum +=  boolc[j]*RadProfile( (b-a)*j*0.25 + a, btt);
	}
	sum += thissum*(b-a)/90.0;
    }

     assert( !std::isnan(sum) && !std::isinf(sum) );

    return sum;
}















