{
    #include "fluxSchemes/HLLFlux.H"

    PtrList<surfaceScalarField> alphaPhiv_hll_m(materials);
    forAll(alphaPhiv_hll_m, mI)
    {
	    alphaPhiv_hll_m.set
	    (
		    mI,
		    new surfaceScalarField
		    (
	           IOobject
               (
	               "alphaPhiv_hll_m"+Foam::name(mI),
                   runTime.timeName(),
                   mesh,
                   IOobject::READ_IF_PRESENT,
                   IOobject::NO_WRITE
               ),
               alphaPhiv_m[mI]
		    )
	    );
    }

    surfaceScalarField phiv_hll = phiv;
    surfaceScalarField phi_hll = phi;
    surfaceVectorField phiUp_hll = phiUp;
    surfaceScalarField phiEp_hll = phiEp;

    #include "fluxSchemes/HLLCFlux.H"

    surfaceScalarField hllc_frac_s = fvc::interpolate(hllc_frac);

    for (int mI=0; mI<materials; mI++) 
    {  
        alphaPhiv_m[mI] = (1.0 - hllc_frac_s) * alphaPhiv_hll_m[mI] + hllc_frac_s * alphaPhiv_m[mI]; 
    }

    phiv  = (1.0 - hllc_frac_s) * phiv_hll  + hllc_frac_s * phiv;
    phi   = (1.0 - hllc_frac_s) * phi_hll   + hllc_frac_s * phi;
    phiUp = (1.0 - hllc_frac_s) * phiUp_hll + hllc_frac_s * phiUp;
    phiEp = (1.0 - hllc_frac_s) * phiEp_hll + hllc_frac_s * phiEp;
}
