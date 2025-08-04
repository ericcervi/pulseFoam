/*-----------------------------------------------------------------------------------------------*\
                        ///        
                        ///                        /////////
  ////////   ///   ///  ///      /////   //////    ///        //////   ///////    //////// /////
  ///   ///  ///   ///  ///    ///      ///  ///   ///       ///  ///       ///   ///  ////  ////
  ///   ///  ///   ///  ///     ///     ////////   ///////   ///  ///   ///////   ///   ///   ///
  ///   ///  ///   ///  ///       ///   ///        ///       ///  ///  ///  ///   ///   ///   ///
  ///   ///  ///   ///  ///        ///  ///        ///       ///  ///  ///  ///   ///   ///   ///
  ////////    ///////    ////  /////     //////    ///        //////    ////////  ///   ///   ///
  ///
  /// 
    
---------------------------------------------------------------------------------------------------
    Built on OpenFOAM v2306
---------------------------------------------------------------------------------------------------
    Copyright (C) 2024 Eric Cervi
---------------------------------------------------------------------------------------------------
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
    pulseFoam

Description
        Density-based multi-material compressible flow solver
        for inertial fusion energy.

\*-----------------------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "pointMesh.C"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolate.H"
#include "fluxSchemes/evalHLL.H"
#include "fluxSchemes/evalHLLC.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes."
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"
    #include "createFields/createCommonFields.H"
    #include "multiMaterialSystem/updatePressure.H"
    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);
    #include "createFields/createCentralFields.H"
    #include "createFields/createHLLCFields.H"

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"

            ++runTime;

            // Do any mesh changes
            mesh.update();
        }

        #include "multiMaterialSystem/updateMultiMatProperties.H"

        // --- Pressure pre-update can improve stability with boiling
        if (boiling)
        {
            #include "multiMaterialSystem/updatePressure.H"
        }

        // --- Directed interpolation of primitive fields onto faces
        #include "updateFields/updateCommonFields.H"

        // Make fluxes relative to mesh-motion
        if (eulRel < 1.0)
        {
            surfaceScalarField meshPhi = fvc::interpolate(W) & mesh.Sf();
            meshPhi.setOriented(false);
            phiv_pos -= meshPhi;
            phiv_neg -= meshPhi;
        }

        c = sqrt(mag(gamma*rPsi));
        c.correctBoundaryConditions();

        if ((fluxScheme == "Kurganov") || (fluxScheme == "Tadmor"))
        {  
            #include "updateFields/updateCentralFields.H"
        }

        if ((fluxScheme == "HLLC") || (fluxScheme == "HLL") || (fluxScheme == "Blended"))
        {  
            #include "updateFields/updateHLLCFields.H"
        }

        // --- ALE mesh routine
        if (eulRel < 1.0)
        {
            //#include "eulEqn.H" 
            forAll(eulerian,cellI)
            {
                eulerian[cellI]=1.0;
            }
            pointInterpolation.interpolate(U, pointU);
    	    pointInterpolation.interpolate(eulerian, pointE);
            pointDisp = pointU*runTime.deltaT();
            pointDisp = pointDisp - eulRel*pointDisp*pointE/max(eulerian);

    	    forAll (pointU, pointI)
    	    {
    	   	    newPoints[pointI] = initialPoints[pointI] + pointDisp[pointI];
    	    }
    		    
    	    mesh.movePoints(newPoints);	
        }

        // --- Update HLLC "star-region" quantities (after mesh motion)
        if ((fluxScheme == "HLLC") || (fluxScheme == "HLL") || (fluxScheme == "Blended"))
        {  
            #include "updateFields/updateHLLCStarFields.H"
        }

        #include "centralCourantNo.H"

        if (LTS)
        {
            #include "setRDeltaT.H"

            ++runTime;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int counter=0; counter<externalCycles; counter++) 
        {            
            // --- Init fluxes       
            #include "createFields/initFluxes.H"

            // --- Calculate fluxes
            if ((fluxScheme == "Kurganov") || (fluxScheme == "Tadmor"))
            {
                #include "fluxSchemes/CentralFlux.H"
            }

            if (fluxScheme == "HLLC")
            {
                #include "fluxSchemes/HLLCFlux.H"
            }

            if (fluxScheme == "HLL")
            {
                #include "fluxSchemes/HLLFlux.H"
            }

            if (fluxScheme == "Blended")
            {
                #include "fluxSchemes/BlendedFlux.H"
            }

            // Make flux for pressure-work absolute
            if (eulRel < 1.0)
            {
                surfaceScalarField meshPhi = fvc::interpolate(W) & mesh.Sf();
                meshPhi.setOriented(false);
                phiEp += meshPhi*(a_pos*p_pos + a_neg*p_neg);
            }

            volScalarField muEff("muEff", rho*nu);
            volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

            // --- Solve material fraction
            #include "continuumMechanics/solvePhaseFractions.H"

            // --- Solve density
            #include "continuumMechanics/solveDensities.H"

            // --- update material properties
            #include "multiMaterialSystem/updateMultiMatProperties.H"

            // --- Solve momentum
            #include "continuumMechanics/solveMomentum.H"

            // --- Solve energy
            #include "continuumMechanics/solveEnergy.H"

            // --- update pressure
            #include "multiMaterialSystem/updatePressure.H"

            // --- update density BC
            //#include "multiMaterialSystem/updateRhoBC.H"
        }

        if (eulRel < 1.0)
        {
    	    initialPoints = newPoints;
        }

        Info << nl << "Mass = " << fvc::domainIntegrate(rho).value() << " kg" << endl;
        Info << nl << "Internal energy = " << fvc::domainIntegrate(rho*Cv*T + 0.5*rho*(U&U)).value() << " J" << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
