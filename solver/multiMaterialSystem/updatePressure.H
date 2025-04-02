{
    psi *= 0.0;
    B0 *= 0.0;
    cavitated *= 0.0;
    boiled *= 0.0;

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

            if (boiling)
            {
                forAll(rho, cellI)
                {
                    if ((p[cellI] > 0.99 * p_cavitation) && (alpha_m[mI][cellI] > 0.01))
                    {
                        if (T[cellI] > 11026.208 / (9.424 - std::log10(std::max(p[cellI] / 133.32, 10.0))))
                        {
                            psi[cellI] += alpha_m[mI][cellI] * ( 1.0 / R_m[mI].value() / T[cellI] - psi_m[mI].value());
                            B0[cellI] -= alpha_m[mI][cellI] * (rho0_m[mI].value() - beta_m[mI].value() * (T[cellI] - T0.value()) - psi_m[mI].value() * p0.value());                    
                            boiled[cellI] = 1.0;
                        }
                    }
                }
            }
        }
    }
    psi.correctBoundaryConditions();
    B0.correctBoundaryConditions();
    boiled.correctBoundaryConditions();

    forAll(B0, cellI)
    {
        if (((rho[cellI] - B0[cellI]) / psi[cellI] < p_cavitation) && cavitating_id < 999)
        {
            cavitated[cellI] = (rho[cellI] 
                               - B0[cellI] 
                               - psi[cellI] * p_cavitation)
                             / (-(rho0_m[cavitating_id].value() - beta_m[cavitating_id].value() * (T[cellI] - T0.value()) - psi_m[cavitating_id].value() * p0.value()) 
                                + 1.0 / R_m[cavitating_id].value() / T[cellI] * p_cavitation
                                - psi_m[cavitating_id].value() * p_cavitation);
                                
            if (cavitated[cellI] > 1.0)
            {
                cavitated[cellI] = 1.0;
            }

            psi[cellI] += cavitated[cellI] * (1.0 / R_m[cavitating_id].value() / T[cellI] - psi_m[cavitating_id].value());

            B0[cellI] -= cavitated[cellI] * (rho0_m[cavitating_id].value() - beta_m[cavitating_id].value() * (T[cellI] - T0.value()) - psi_m[cavitating_id].value() * p0.value());
        }
    }
    psi.correctBoundaryConditions();
    B0.correctBoundaryConditions();
    cavitated.correctBoundaryConditions();

    p = (rho - B0)/psi;
    p.correctBoundaryConditions();
}
