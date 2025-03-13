{
    psi *= 0.0;
    B0 *= 0.0;
    cavitated *= 0.0;
    for (int mI=0; mI<materials; mI++) 
    {
        if(EOS_m[mI]==0)
        {
            psi += alpha_m[mI] / R_m[mI] / T;
        }
        if(EOS_m[mI]==1)
        {
            psi += alpha_m[mI] * psi_m[mI];
            B0 += alpha_m[mI] * (rho0_m[mI] - beta_m[mI] * (T - T0) - psi_m[mI] * p0);
        }
    }

    forAll(B0, cellI)
    {
        if ((rho[cellI] - B0[cellI]) / psi[cellI] < p_cavitation)
        {
            psi[cellI] *= 0.0;
            for (int mI=0; mI<materials; mI++) 
            {
                psi[cellI] += alpha_m[mI][cellI] / R_m[mI].value() / T[cellI];
            }
            if (psi[cellI] <= 0.0)
            {
                psi[cellI] = 1.0e-9;
            }
            B0[cellI] = 0.0;
            cavitated[cellI] = 1.0;
        }
    }
    psi.correctBoundaryConditions();
    B0.correctBoundaryConditions();
    cavitated.correctBoundaryConditions();

    p = (rho - B0)/psi;
    p.correctBoundaryConditions();
}
