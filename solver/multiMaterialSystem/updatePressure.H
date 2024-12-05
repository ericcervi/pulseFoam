{
    psi *= 0.0;
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
            B0 += alpha_m[mI] * (rho0_m[mI] - beta_m[mI]*(T - T0) - psi_m[mI]*p0);
        }
    }

    volScalarField diff = B0 * 0.0;
    forAll(B0, cellI)
    {
        if (psi[cellI] <= 0.0)
        {
            psi[cellI] = 1e-9;
        }
        if (B0[cellI] >= rho[cellI])
        {
            diff[cellI] = B0[cellI] - rho[cellI] + 1e-9;
        }
    }
    diff.correctBoundaryConditions();
    psi.correctBoundaryConditions();
    B0.correctBoundaryConditions();

    p = (rho - B0 + diff)/psi;
    p.correctBoundaryConditions();
}
