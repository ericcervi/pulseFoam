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
//////////// TEST 
            #include "multiMaterialSystem/updateMultiMatProperties.H"
 
/*psi *= 0.0;
B0 *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    if(EOS_m[mI]==0)
    {
        psi += alpha_m[mI] / R_m[mI] / T;
    }
    if(EOS_m[mI]==1)
    {
        psi += alpha_m[mI] * psi_m[mI];
        B0 +=  alpha_m[mI] * (rho0_m[mI] - beta_m[mI]*(T - T0) - psi_m[mI]*p0);
    }
}

forAll(psi, cellI)
{
    if (psi[cellI] <= 0.0)
    {
        psi[cellI] = 1e-9;
    }
}
psi.correctBoundaryConditions();
B0.correctBoundaryConditions();*/           


///////////////
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

        if (fluxScheme == "HLLC")
        {  
            #include "updateFields/updateHLLCFields.H"
        }

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
        if (fluxScheme == "HLLC")
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

            // Make flux for pressure-work absolute
            if (eulRel < 1.0)
            {
                surfaceScalarField meshPhi = fvc::interpolate(W) & mesh.Sf();
                meshPhi.setOriented(false);
                phiEp += meshPhi*(a_pos*p_pos + a_neg*p_neg);
            }

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
            #include "multiMaterialSystem/updateRhoBC.H"
        }

        if (eulRel < 1.0)
        {
    	    initialPoints = newPoints;
        }

        Info << nl << "Mass = " << fvc::domainIntegrate(rho).value() << " kg" << endl;
        Info << nl << "Enthalpy = " << fvc::domainIntegrate(rho*Cv*T + p).value() << " kg" << endl;
        Info << nl << "Total Enthalpy = " << fvc::domainIntegrate(rho*Cv*T + p + 0.5*rho*(U&U)).value() << " kg" << endl;
        Info << nl << "Internal energy = " << fvc::domainIntegrate(rho*Cv*T + 0.5*rho*(U&U)).value() << " kg" << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
