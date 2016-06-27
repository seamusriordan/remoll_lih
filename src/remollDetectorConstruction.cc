

#include "remollDetectorConstruction.hh"
#include "remollGenericDetector.hh"
#include "remollBeamTarget.hh"
#include "remollGlobalField.hh"
#include "remollRun.hh"
#include "remollRunData.hh"
#include "remollIO.hh"

#include "TGeoManager.h"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "globals.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

// GDML export
#include "G4GDMLParser.hh"

//visual
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#define __DET_STRLEN 200
#define __MAX_DETS 5000

remollDetectorConstruction::remollDetectorConstruction() {
    // Default geometry file
    fDetFileName = "geometry_sculpt/mollerMother.gdml";


    CreateGlobalMagneticField();
    fIO = NULL;
    fGDMLParser = NULL;
}

remollDetectorConstruction::~remollDetectorConstruction() {
}

G4VPhysicalVolume* remollDetectorConstruction::Construct() {

    double z,a,density;
    int nelements;

    G4Element* Al = new G4Element("Aluminum", "Al", z=13 , a=27*g/mole);
    G4Material* Alu_mat = new G4Material("Alu_Mat", 2.7*g/cm3, nelements=1);
    Alu_mat->AddElement(Al, 1);


    G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

    G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
    Air->AddElement(N, 70.*perCent);
    Air->AddElement(O, 30.*perCent);

    double world_x, world_y, world_z;
    world_x = world_y = world_z = 275*cm;
    G4Box* world_box = new G4Box("World",world_x,world_y,world_z);
    G4LogicalVolume* world_log
        = new G4LogicalVolume(world_box,Air,"World",0,0,0);

    double target_thick = 0.2*mm;


    G4Tubs* target_center= new G4Tubs("target_center", 0, 1.25*cm, target_thick/2, 0.0, 360*deg);
    G4LogicalVolume *target_center_log = new G4LogicalVolume(target_center, Alu_mat, "target_center_log", 0, 0, 0);
    G4VPhysicalVolume* target_center_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),target_center_log,"target_center_phys",world_log,false,0);


    G4Tubs* target_upwall= new G4Tubs("target_upwall", 0, 1.25*cm, 100*um/2, 0.0, 360*deg);
    G4LogicalVolume *target_upwall_log = new G4LogicalVolume(target_upwall, Alu_mat, "target_upwall_log", 0, 0, 0);
    G4VPhysicalVolume* target_upwall_phys = new G4PVPlacement(0,G4ThreeVector(0,0,target_thick/2+50*um),target_upwall_log,"target_upwall_phys",world_log,false,0);

    G4Tubs* target_downwall= new G4Tubs("target_downwall", 0, 1.25*cm, 100*um/2, 0.0, 360*deg);
    G4LogicalVolume *target_downwall_log = new G4LogicalVolume(target_downwall, Alu_mat, "target_downwall_log", 0, 0, 0);
    G4VPhysicalVolume* target_downwall_phys = new G4PVPlacement(0,G4ThreeVector(0,0,target_thick/2+50*um),target_downwall_log,"target_downwall_phys",world_log,false,0);

     G4VPhysicalVolume* worldVolume
         = new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);


     // Set target volume for generator

     remollBeamTarget *beamtarg = remollBeamTarget::GetBeamTarget();
     beamtarg->Reset();
     // Choose where these come from
     beamtarg->AddVolume(target_center_phys);
//     beamtarg->AddVolume(target_upwall_phys);
//     beamtarg->AddVolume(target_downwall_phys);

    G4Box* detector = new G4Box("detector", 1*m, 1*m, 1*mm );
    G4LogicalVolume *detector_log = new G4LogicalVolume(detector, Alu_mat, "target_upwall_log", 0, 0, 0);
    G4VPhysicalVolume* detector_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),detector_log,"detector_phys",world_log,false,0);

     G4SDManager* SDman = G4SDManager::GetSDMpointer();
     SDman->AddNewDetector(new remollGenericDetector("mydet", 1));


    return worldVolume;
}

G4int remollDetectorConstruction::UpdateCopyNo(G4VPhysicalVolume* aVolume,G4int index){  

  //if (aVolume->GetLogicalVolume()->GetNoDaughters()==0 ){
      aVolume->SetCopyNo(index);
      index++;
      //}else {
    for(int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++){
      index = UpdateCopyNo(aVolume->GetLogicalVolume()->GetDaughter(i),index);
    }
    //}

  return index;
};

void remollDetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth)
{
  for(int isp=0;isp<depth;isp++)
  { G4cout << "  "; }
  //aVolume->SetCopyNo(1);
  G4cout << aVolume->GetName() << "[" << aVolume->GetCopyNo() << "] "
         << aVolume->GetLogicalVolume()->GetName() << " "
         << aVolume->GetLogicalVolume()->GetNoDaughters() << " "
         << aVolume->GetLogicalVolume()->GetMaterial()->GetName() << " "
	 << G4BestUnit(aVolume->GetLogicalVolume()->GetMass(true),"Mass");
  if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
  {
    G4cout << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()
                            ->GetFullPathName();
  }
  G4cout << G4endl;
  for(int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++)
  { DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),depth+1); }
}

void remollDetectorConstruction::CreateGlobalMagneticField() {
    fGlobalField = new remollGlobalField();

    fGlobalFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fGlobalFieldManager->SetDetectorField(fGlobalField);
    fGlobalFieldManager->CreateChordFinder(fGlobalField);

    return;
} 

void remollDetectorConstruction::SetDetectorGeomFile(const G4String &str){
    fDetFileName = str;
}
