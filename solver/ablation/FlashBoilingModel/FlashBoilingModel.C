/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    initializeAblation

Group
    grpPreProcessingUtilities

Description
    Initialize simulation in order to model the ablation
    of a liquid wall caused by the energy deposition of a
    radiation source.
    One needs to set the position of the radiation source,
    the name of the patch representing the wall, and some 
    parameters characterizing the wall and the radiation 
    source.


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "patchExprDriver.H"
#include "timeSelector.H"
// #include "readFields.H"
#include "FaceShading.H"
#include "HeatDepositionTable.H"
#include "FlashBoilingModel.H"

namespace Foam
{
    defineTypeNameAndDebug(FlashBoilingModel, 0);
}




Foam::FlashBoilingModel::FlashBoilingModel(const Time& runTime, const fvMesh& mesh):
mesh_(mesh),
initDict_
(
    IOobject
    (
        "initAblationDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
origin_(initDict_.get<Vector<scalar>>("origin")),
found_(false),
pPtr_(nullptr),
TPtr_(nullptr),
rhoPtr_(nullptr),
faceZoneIDs_(0),
cellZoneIDs_(0),
faceShadingPtr_(nullptr)
{
    Info << "Flash boiling model constructor " << endl;

    // Setting a list of field names whose existence to be checked
    wordList fieldNames(initDict_.get<wordList>("fields"));
    

    for(auto& fieldNameI : fieldNames)
    {
        // dummy IOobject to check whether the field named as the entry in "field" exists
        IOobject fieldHeaderI
        (
            fieldNameI,   // name of the IOobject
            mesh.thisDb().time().timeName(), // name of the file
            mesh.thisDb(), // object registry
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );    
        
        const bool headiOk = fieldHeaderI.typeHeaderOk<IOdictionary>(false);

        if (!headiOk)
        {

            WarningInFunction
                        << "Requested field to change " << fieldNameI
                        << " does not exist in " << fieldHeaderI.path() << endl;

            found_ = false;

            // Do something for error control
        }
        else 
        {
            Info << fieldNameI << " found." << endl;
            found_ = true;

        }
    }

    if (found_)
    {
        Info<< "Reading field p\n" << endl;
        pPtr_ = std::make_unique<volScalarField>
        (
            IOobject
            (
                fieldNames[0],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Reading field T\n" << endl;
        TPtr_ = std::make_unique<volScalarField>
        (
            IOobject
            (
                fieldNames[1],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );


        Info<< "Reading field rho\n" << endl;
        rhoPtr_ = std::make_unique<volScalarField>
        (
            IOobject
            (
                fieldNames[2],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        
        readZones();

        faceShadingPtr_ = std::make_unique<FaceShading>(mesh, labelList(1,-1) , faceZoneIDs_, origin_);
    }


}

bool Foam::FlashBoilingModel::readZones()
{
    // faceZoneMesh is a list of faceZones
    // The indexes are global unlike for patches
    const faceZoneMesh& fzm = mesh_.faceZones();
    
    // Find faceZones names
    const wordList faceZoneNames(initDict_.get<wordList>("faceZones"));
    
    // Find zoneIDs associated with patchNames
    for(const auto& namei : faceZoneNames)
    {
        faceZoneIDs_.append(fzm.findZoneID(namei));
    }

    // cellZoneMesh is a list of cellZones
    const cellZoneMesh& czm = mesh_.cellZones();
    
    // Find faceZones names
    const wordList cellZoneNames(initDict_.get<wordList>("cellZones"));
    
    // Find zoneIDs associated with patchNames
    for(const auto& namei : cellZoneNames)
    {
        cellZoneIDs_.append(czm.findZoneID(namei));
    }

    return true;
}

void Foam::FlashBoilingModel::ablate()
{
    //
    // Model parameters
    //
    const scalar Tcrit = initDict_.get<scalar>("Tcrit");
    
    // Explosive boiling threshold
    const scalar qExp = initDict_.get<scalar>("qExp");
    
    // Average volumetric energy deposition
    const scalar qAv = initDict_.get<scalar>("qAv");
        
    // specific heat of the gas
    const scalar cv = initDict_.get<scalar>("cv");
    
    // Liquid density
    const scalar rhol = initDict_.get<scalar>("rhol");
    
    // Molar mass
    const scalar MM = initDict_.get<scalar>("MM");
    
    // Specific gas constant
    const scalar Rsp = 8314/MM;
    
    // Name of the subdict associated with the heat deposition table
    const word tableName("heatDeposition");

    //
    // Get references to fields
    // 

    volScalarField& p   = *pPtr_;
    volScalarField& T   = *TPtr_;
    volScalarField& rho = *rhoPtr_;

    // Only for debugging purposes
    label  nFaces = 0;
    scalar ablatedMass = 0.0;
    scalar totalArea = 0.0;
    // scalar totalF = 0.0;
    scalar averageAblatedDepth = 0.0;


    //
    // Topology info
    //

    const labelList& hitFaces(faceShadingPtr_->visibleFaces());

    const cellZoneMesh& czm = mesh_.cellZones();
    const faceZoneMesh& fzm = mesh_.faceZones();


    // Currently working on one faceZone only
    const labelList& hitFacesFront = fzm[0].frontCells();
    const labelList& hitFacesBack = fzm[0].backCells();
    
    Info << "hitFace[0] " << hitFaces[0] << endl;
    
    for (int i = 0; i < hitFaces.size() ; ++i)
    {  
        // whichPatchFace gives back a pair of labels representing the patchID
        // and the local face index on the patch.
        /*label patchI  = pbm.whichPatchFace(hitFaces[i])[0];
        
        // I need to use local indexing
        label faceI = hitFaces[i] - pbm[patchI].start();*/


        label faceI = hitFaces[i];
        label zoneI = fzm.whichZone(faceI);

        // map from global face index to local 
        const Map<label>& labelMap = fzm[zoneI].lookupMap();

        label localFaceI = labelMap[faceI];

        label gasCellI = -1; // Initialize to invalid label

        // Check if the neighbouring gas is front or back using the provided
        // cellZone name. cellZoneIDs[0] must be the name of the gas region
        if (czm.whichZone(hitFacesFront[localFaceI]) == cellZoneIDs_[0])
        {
            gasCellI = hitFacesFront[localFaceI];
        }
        else
        {
            gasCellI = hitFacesBack[localFaceI];
        }

        
        const primitiveFacePatch& pp = fzm[zoneI]();
        const scalar& magSfI = pp.magFaceAreas()[localFaceI];
        const vector& normalI = pp.faceNormals()[localFaceI];
        const vector& zoneFaceCentreI = pp.faceCentres()[localFaceI];

        const vector sourceToZoneFaceI = zoneFaceCentreI - origin_;
        const vector rayUnitVectorI = sourceToZoneFaceI/mag(sourceToZoneFaceI);
        const scalar magSourceToZoneFaceI = mag(sourceToZoneFaceI);

        scalar viewFactorI = mag(normalI & rayUnitVectorI) * 0.5 * 0.5 / (magSourceToZoneFaceI * magSourceToZoneFaceI); // abs() acts wierd
            
        //scalar scalarProduct = normalI & rayUnitVectorI;
        
        // Check boiling condition
        if (qAv * viewFactorI > qExp)
        {
            Function1Types::HeatDepositionTable tableI(tableName, initDict_, qAv * viewFactorI);
            scalar ablationDepthI = tableI.reverseValue(qExp);
            
            scalar ablatedMassI = rhol * ablationDepthI * 1e-6 * magSfI;
            
            // Spread the estimated ablated mass over the owner cell volume            
            rho[gasCellI] = ablatedMassI/mesh_.V()[gasCellI];
            

            scalar depositedVolHeatI = tableI.integrate(0 , ablationDepthI) * magSfI * 1e-6;
            depositedVolHeatI /= mesh_.V()[gasCellI];


            T[gasCellI] = 0.844 * Tcrit + depositedVolHeatI/rho[gasCellI]/cv;
            
            p[gasCellI] = rho[gasCellI] * Rsp * T[gasCellI];

            ablatedMass += ablatedMassI;
            totalArea += magSfI;
            averageAblatedDepth += ablationDepthI;
            nFaces++;
        }
         
    }

    Info << "Number of ablated faces " << nFaces << endl;
    Info << "Estimated ablatedMass " << ablatedMass << " kg"<< endl;
    Info << "average ablatedDepth " << averageAblatedDepth/nFaces << " micron"<< endl;
    Info << "Total area " << totalArea << " m2" << endl;
    
    // // Write fields to disk
    T.write();
    rho.write();
    p.write();
}

// Memory cleanup handled by smart pointers
Foam::FlashBoilingModel::~FlashBoilingModel()
{}



// ************************************************************************* //
