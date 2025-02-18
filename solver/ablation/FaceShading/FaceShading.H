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
    Foam::FaceShading

Description
    Helper class to calculate visible faces for point source to face illumination.
    Similar to FaceShading in OpenFOAM except the source doesn't come from a unique 
    direction

SourceFiles
    FaceShading.C

\*---------------------------------------------------------------------------*/

#ifndef faceTracing_H
#define faceTracing_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "triSurface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FaceShading
{
// Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Patches to check for visibility
        const labelList patchIDs_;


        labelList zoneIDs_;

        //- Direction
        point sourcePoint_;

        //- Faces directly hit by vector direction
        labelList visibleFaces_;


        //- Faces directed towards the source
        labelList directHitFaces_;


    // Private Member Functions

        //- Calculate visible faces from patches
        void calculateFromPatches();

        //- Calculate visible faces from zones
        void calculateFromZones();

        //- No copy construct
        FaceShading(const FaceShading&) = delete;

        //- No copy assignment
        void operator=(const FaceShading&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("FaceShading");


    // Constructors

        //- Construct from mesh and vector to 'sun'. All uncoupled patches
        //  are checked for visibility. faceZones are ignored.
        FaceShading
        (
            const fvMesh& mesh,
            const label& patchID,
            const point& sourcePoint
        );


        FaceShading
        (
            const fvMesh& mesh,
            const labelList& patchIDs,
            const labelList& zoneIDs,
            const point& sourcePoint
        );
        

    


    //- Destructor
    ~FaceShading() = default;


    // Member Functions

        void selectFaces 
        (
            bool useNormal,
            const label& patchID,
            labelList& faceIDs,
            bitSet& flipMap
        ) const;

        // Triangulate patch
        triSurface triangulate
        (
            const labelUList& faceIDs,
            const bitSet& flipMap
        ) const;

        // Triangulate patch
        triSurface triangulateZones
        (
            const labelUList& faceIDs,
            const bitSet& flipMap
        ) const;


        void selectFaces 
        (
            bool useNormal,
            const labelList& patchIDs,
            labelList& faceIDs,
            bitSet& flipMap
        ) const;


        void selectFacesFromZones 
        (
            bool useNormal,
            const labelList& zoneIDs,
            labelList& faceIDs,
            bitSet& flipMap
        ) const;

        // Access

        //- Const access to rayStartFaces
        const labelList& visibleFaces() const
        {
            return visibleFaces_;
        }

        // Just for debugging
        const labelList& directHitFaces() const
        {
            return directHitFaces_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
