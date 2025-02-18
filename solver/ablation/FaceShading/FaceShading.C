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

\*---------------------------------------------------------------------------*/

#include "FaceShading.H"
#include "distributedTriSurfaceMesh.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FaceShading, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::FaceShading::calculateFromPatches()
{
    labelList hitFacesIds;
    bitSet hitFacesFlips;
    selectFaces
    (
        true, // use normal to do first filtering
        patchIDs_,
        hitFacesIds,
        hitFacesFlips
    );

    Info<< "Number of 'potential' direct hits : "
        << returnReduce(hitFacesIds.size(), sumOp<label>()) << endl;

    // Find potential obstructions. Include all faces that might potentially
    // block (so ignore normal by setting useNormal param to false)
    labelList blockingFacesIds;
    bitSet blockingFacesFlips;
    selectFaces
    (
        false,   
        patchIDs_,
        blockingFacesIds,
        blockingFacesFlips
    );

    if (debug)
    {
        Info << "Number of blocking faces " << blockingFacesIds.size() << endl; 
    }

    // Triangulate faces to use the topology search algorithm of OpenFOAM
    const triSurface localSurface = triangulate
    (
        blockingFacesIds,
        blockingFacesFlips
    );

    // Create distributedTriSurfaceMesh
    Random rndGen(653213);

    // Determine mesh bounding boxes:
    List<treeBoundBox> meshBb
    (
        1,
        treeBoundBox(mesh_.points()).extend(rndGen, 1e-3)
    );

    // // Dummy bounds dictionary
    dictionary dict;
    dict.add("bounds", meshBb);
    dict.add
    (
        "distributionType",
        distributedTriSurfaceMesh::distributionTypeNames_
        [
            distributedTriSurfaceMesh::FROZEN
        ]
    );
    dict.add("mergeDistance", SMALL);

    // surfacesMesh represents the triangulated surface mesh composed of 
    // all the possible blocking faces identified in the mesh.
    // It is then used to find the (possible) intersection between rays
    // and faces.
    distributedTriSurfaceMesh surfacesMesh
    (
        IOobject
        (
            "opaqueSurface.stl",
            mesh_.time().constant(),    // directory
            "triSurface",               // instance
            mesh_.time(),               // registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        localSurface,
        dict
    );
    
    // Shoot rays
    {
        DynamicField<point> start(hitFacesIds.size());
        DynamicList<label> endIndex(start.size());
        DynamicField<point> end(hitFacesIds.size());

        const pointField& faceCentres = mesh_.faceCentres();

        forAll(hitFacesIds, i)
        {
            const label facei = hitFacesIds[i];
            const point& fc = faceCentres[facei];
            //Info << "fc " << faceCentres[facei] << endl;
            vector sourceToPatchFacei = fc - sourcePoint_;

            // Set the end point of the ray just before the face
            // to not have a numerical self-shading of a face that
            // is actually visible to the source.
            // The aim of the code is to detect obstructions from 
            // OTHER faces.
            end.append(fc - 1e-4 * sourceToPatchFacei);
            endIndex.append(facei);
            start.append(sourcePoint_);
        }

        List<pointIndexHit> hitInfo(endIndex.size());

        // Update hits
        surfacesMesh.findLine(start, end, hitInfo);

        label nVisible = 0;

        // loop over rays and check visibility
        // between a ray and a face.
        // If a hit is found, then the face is not visible
        forAll(hitInfo, rayI)
        {
           if (!hitInfo[rayI].hit())
            {
                nVisible++;
            }
            
        }
        
        visibleFaces_.setSize(nVisible);
        nVisible = 0;

        forAll(hitInfo, rayI)
        {
           if (!hitInfo[rayI].hit())
            {
                visibleFaces_[nVisible++] = endIndex[rayI];
            }   
        }

        start.clear();
        endIndex.clear();
        end.clear();
        
    }

    Info << "Number of visible faces " << visibleFaces_.size() << endl;
}


void Foam::FaceShading::calculateFromZones()
{
    labelList hitFacesIds;
    bitSet hitFacesFlips;

    // Do not apply first filtering on zones
    // based on their orientation because internal faces
    // can have arbitrary orientations.
    selectFacesFromZones
    (
        false, 
        zoneIDs_,
        hitFacesIds,
        hitFacesFlips
    );

    Info<< "Number of 'potential' direct hits : "
        << returnReduce(hitFacesIds.size(), sumOp<label>()) << endl;

    // Find potential obstructions. Include all faces that might potentially
    // block (so ignore normal by setting useNormal param to false)
    labelList blockingFacesIds;
    bitSet blockingFacesFlips;
    selectFacesFromZones
    (
        false,   
        zoneIDs_,
        blockingFacesIds,
        blockingFacesFlips
    );

    Info << "Number of blocking faces : " << blockingFacesIds.size() << endl; 
    
    const triSurface localSurface = triangulateZones
    (
        blockingFacesIds,
        blockingFacesFlips
    );

    // Create distributedTriSurfaceMesh
    Random rndGen(653213);

    // Determine mesh bounding boxes:
    List<treeBoundBox> meshBb
    (
        1,
        treeBoundBox(mesh_.points()).extend(rndGen, 1e-3)
    );

    // // Dummy bounds dictionary
    dictionary dict;
    dict.add("bounds", meshBb);
    dict.add
    (
        "distributionType",
        distributedTriSurfaceMesh::distributionTypeNames_
        [
            distributedTriSurfaceMesh::FROZEN
        ]
    );
    dict.add("mergeDistance", SMALL);

    // surfacesMesh represents the triangulated surface mesh composed of 
    // all the possible blocking faces identified in the mesh.
    // It is then used to find the (possible) intersection between rays
    // and faces.
    distributedTriSurfaceMesh surfacesMesh
    (
        IOobject
        (
            "opaqueSurface.stl",
            mesh_.time().constant(),    // directory
            "triSurface",               // instance
            mesh_.time(),               // registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        localSurface,
        dict
    );
    
    // Shoot rays
    {
        DynamicField<point> start(hitFacesIds.size());
        DynamicList<label> endIndex(start.size());
        DynamicField<point> end(hitFacesIds.size());

        const pointField& faceCentres = mesh_.faceCentres();

        forAll(hitFacesIds, i)
        {
            const label facei = hitFacesIds[i];
            const point& fc = faceCentres[facei];
            vector sourceToPatchFacei = fc - sourcePoint_;

            end.append(fc - 1e-4 * sourceToPatchFacei); // to be checked for multiple faceZones
            endIndex.append(facei);
            start.append(sourcePoint_);
        }

        List<pointIndexHit> hitInfo(endIndex.size());
        surfacesMesh.findLine(start, end, hitInfo);

        //Info << "hitInfo " << hitInfo << endl;

        label nVisible = 0;
        forAll(hitInfo, rayI)
        {
           if (!hitInfo[rayI].hit())
            {
                nVisible++;
            }
            
        }
        
        visibleFaces_.setSize(nVisible);
        nVisible = 0;

        forAll(hitInfo, rayI)
        {
           if (!hitInfo[rayI].hit())
            {
                visibleFaces_[nVisible++] = endIndex[rayI];
            }   
        }

        start.clear();
        endIndex.clear();
        end.clear();
        
    }

    Info << "Number of visible faces " << visibleFaces_.size() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FaceShading::FaceShading
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& zoneIDs,
    const point& sourcePoint 
)
:
    mesh_(mesh),
    patchIDs_(patchIDs),
    zoneIDs_(zoneIDs),
    sourcePoint_(sourcePoint),
    visibleFaces_(0),
    directHitFaces_(0)
{
    calculateFromZones();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FaceShading::selectFaces
(
    const bool useNormal,
    const labelList& patchIDs,
    labelList& faceIDs,
    bitSet& flipMap
) const
{
    // reference to boundary mesh 
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    bitSet isSelected(mesh_.nFaces());
    DynamicList<label> dynFaces(mesh_.nBoundaryFaces());
    bitSet isFaceFlipped(mesh_.nFaces());


    // Loop over the faces of the given patch 

    for(const label patchi : patchIDs)
    {
        const polyPatch pp = pbm[patchi];

        const vectorField& n = pp.faceNormals();

        forAll(n, i)
        {
            const label facei = i + pp.start();
            const vector sourceToPatchFacei(mesh_.faceCentres()[facei] - sourcePoint_);

            if (!useNormal || (sourceToPatchFacei & n[i] ) > 0)
            {
                isSelected.set(facei);
                isFaceFlipped[facei] = false;
                dynFaces.append(facei);
            }
        } 
    }

    faceIDs = std::move(dynFaces);
    flipMap = bitSet(isFaceFlipped, faceIDs);
}

// There's a problem with the orientation of face zones
void Foam::FaceShading::selectFacesFromZones
(
    const bool useNormal,
    const labelList& zoneIDs,
    labelList& faceIDs,
    bitSet& flipMap
) const
{
    // Select from faceZones
    const auto& fzs = mesh_.faceZones();


    bitSet isSelected(mesh_.nFaces());
    DynamicList<label> dynFaces(mesh_.nBoundaryFaces());
    bitSet isFaceFlipped(mesh_.nFaces());

    for (const label zonei : zoneIDs)
    {
        const auto& fz = fzs[zonei];
        const primitiveFacePatch& pp = fz();
        // I need to make sure the primitive patch is oriented how I want
        // I need to define a new vector field pointing towards the gas 
        const vectorField& n = pp.faceNormals();


        forAll(n, i)
        {
            const label facei = fz[i] ;
            const vector sourceToPatchFacei(mesh_.faceCentres()[facei] - sourcePoint_);

            scalar orientation = sourceToPatchFacei & n[i] ;


            if(!useNormal /*|| (sourceToPatchFacei & n[i] ) < 0*/)
            {
                isSelected.set(facei);
                dynFaces.append(facei);
                isFaceFlipped[facei] = true;
            }
        }
    }

    faceIDs = std::move(dynFaces);
    flipMap = bitSet(isFaceFlipped, faceIDs);
}

Foam::triSurface Foam::FaceShading::triangulate
(
    const labelUList& faceIDs,
    const bitSet& flipMap
) const
{
    if (faceIDs.size() != flipMap.size())
    {
        FatalErrorInFunction << "Size problem :"
            << "faceIDs:" << faceIDs.size()
            << "flipMap:" << flipMap.size()
            << exit(FatalError);
    }

    const auto& points = mesh_.points();
    const auto& faces = mesh_.faces();
    const auto& bMesh = mesh_.boundaryMesh();

    geometricSurfacePatchList surfPatches(bMesh.nNonProcessor());
    labelList patchID(mesh_.nFaces(), -1);
    {
        label newPatchi = 0;
        for (label patchi = 0; patchi < bMesh.nNonProcessor(); ++patchi)
        {
            const auto& pp = bMesh[patchi];

            surfPatches[newPatchi] = geometricSurfacePatch
            (
                pp.name(),
                newPatchi,
                pp.type()
            );
            SubList<label>
            (
                patchID,
                pp.size(),
                pp.start()
            ) = newPatchi;

            newPatchi++;
        }
    }


    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(2*faceIDs.size());

    // Work array
    faceList triFaces;

    forAll(faceIDs, i)
    {
        const label facei = faceIDs[i];
        const bool flip = flipMap[i];
        const label patchi = patchID[facei];
        const face& f = faces[facei];

        // Triangulate face
        triFaces.setSize(f.nTriangles(points));
        label nTri = 0;
        f.triangles(points, nTri, triFaces);

        for (const face& f : triFaces)
        {
            if (!flip)
            {
                triangles.append(labelledTri(f[0], f[1], f[2], patchi));
            }
            else
            {
                triangles.append(labelledTri(f[0], f[2], f[1], patchi));
            }
        }
    }

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh_.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().transfer(surfPatches);

    return surface;
}


Foam::triSurface Foam::FaceShading::triangulateZones
(
    const labelUList& faceIDs,
    const bitSet& flipMap
) const
{
    if (faceIDs.size() != flipMap.size())
    {
        FatalErrorInFunction << "Size problem :"
            << "faceIDs:" << faceIDs.size()
            << "flipMap:" << flipMap.size()
            << exit(FatalError);
    }

    const auto& points = mesh_.points();
    const auto& faces = mesh_.faces();
    const auto& bMesh = mesh_.boundaryMesh();
    const auto& fzs = mesh_.faceZones();


    // Patching of surface:
    // - non-processor patches
    // - faceZones
    // Note: check for faceZones on boundary? Who has priority?
    geometricSurfacePatchList surfPatches(fzs.size());

    // Debug
    //Info << "surfPatches size " << surfPatches.size() << endl;
    // End debug

    labelList patchID(mesh_.nFaces(), -1);
    {
        label newPatchi = 0;
        for (const auto& fz : fzs)
        {
            surfPatches[newPatchi] = geometricSurfacePatch
            (
                fz.name(),
                newPatchi,
                fz.type()
            );
            UIndirectList<label>(patchID, fz) = newPatchi;

            newPatchi++;
        }
    }

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(2*faceIDs.size());

    // Work array
    faceList triFaces;

    // Debug
    //Info << "faceIDs " << faceIDs << endl;
    // End debug

    forAll(faceIDs, i)
    {
        const label facei = faceIDs[i];
        const bool flip = flipMap[i];
        const label patchi = patchID[facei];
        const face& f = faces[facei];

        // Triangulate face
        triFaces.setSize(f.nTriangles(points));
        label nTri = 0;
        f.triangles(points, nTri, triFaces);

        for (const face& f : triFaces)
        {
            if (!flip)
            {
                triangles.append(labelledTri(f[0], f[1], f[2], patchi));
            }
            else
            {
                triangles.append(labelledTri(f[0], f[2], f[1], patchi));
            }
        }
    }

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh_.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    //Info << "surface " << surface << endl;

    // Add patch names to surface
    surface.patches().transfer(surfPatches);

    return surface;
}

// ************************************************************************* //
