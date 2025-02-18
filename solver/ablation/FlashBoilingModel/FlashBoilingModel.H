/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

Class
    Foam::FlashBoilingModel

Description
    Ablation modeling based on flash boiling

SourceFiles
    FlashBoilingModel.C

\*---------------------------------------------------------------------------*/

#ifndef flashBoiling_H
#define flashBoiling_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// Just for testing

class FlashBoilingModel
{
    // Private data

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Dictionary
    IOdictionary initDict_;

    //- Radiation source location
    Vector<scalar> origin_;

    //- Check if fields are found
    bool found_;

    //- Pressure
    std::unique_ptr<volScalarField> pPtr_;

    //- Temperature
    std::unique_ptr<volScalarField> TPtr_;

    //- Gas density
    std::unique_ptr<volScalarField> rhoPtr_;

    //- faceZones IDs
    labelList faceZoneIDs_;

    //- cellZones IDs
    labelList cellZoneIDs_;

    //- pointer to ray-tracing routine
    std::unique_ptr<FaceShading> faceShadingPtr_;


    // Private functions
    bool readZones();


public:
    
    // Declare name of the class and its debug switch
    ClassName("FaceShading");


    FlashBoilingModel(const Time& runTime, const fvMesh& mesh);
    ~FlashBoilingModel();

    void ablate();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //