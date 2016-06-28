#ifndef __REMOLLGENMOTT_HH 
#define __REMOLLGENMOTT_HH 
/*!
 * Mott elastic event generator
 *
 * Seamus Riordan
 * June 27, 2016
 *
 * Based heavily on previous work from mollersim
*/

#include "remollVEventGen.hh"

class remollBeamTarget;

class remollGenMott : public remollVEventGen {
    public:
	remollGenMott();
	~remollGenMott();

    private:
	void SamplePhysics(remollVertex *, remollEvent *);

	G4double RadProfile(G4double,G4double);
	G4double EnergNumInt(G4double,G4double,G4double);

	remollBeamTarget *fBeamTarg;
};

#endif//__REMOLLGENPELASTIC_HH 
