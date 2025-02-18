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

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
int main(int argc, char *argv[])
{
    // No -constant, no special treatment for 0/
    timeSelector::addOptions(false);

    // "setRootCase.H" defines the argument list object
    #include "setRootCase.H"

    // Creates the time object, self-explanatory
    #include "createTime.H"

    #include "createNamedMesh.H"

    Info << "\nSetting IO object that reads dictionary\n" << endl;

    // dictionary name to look at. Name of the file
    const word dictName("initAblationDict");

    FlashBoilingModel flashBoiling(runTime, mesh);

    flashBoiling.ablate();
    
    Info << "\n End " << endl;

    return 0;
}


// ************************************************************************* //
