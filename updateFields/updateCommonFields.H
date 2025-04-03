{
    for (int mI=0; mI<materials; mI++) 
    {
        if (materials > 1)
        {
            alpha_pos_m[mI] = interpolate(alpha_m[mI], pos, rho.name());
            alpha_neg_m[mI] = interpolate(alpha_m[mI], neg, rho.name());
        }

        rho_pos_m[mI] = interpolate(rho_m[mI], pos, rho.name());
        rho_neg_m[mI] = interpolate(rho_m[mI], neg, rho.name());
    }

    rho_pos = interpolate(rho, pos);
    rho_neg = interpolate(rho, neg);
        
    rhoU_pos = interpolate(rhoU, pos, U.name());
    rhoU_neg = interpolate(rhoU, neg, U.name());

    rPsi = 1.0/psi;
    rPsi.correctBoundaryConditions();
    rPsi_pos = interpolate(rPsi, pos, T.name());
    rPsi_neg = interpolate(rPsi, neg, T.name());

    B0_pos = interpolate(B0, pos, rho.name());
    B0_neg = interpolate(B0, neg, rho.name());

    e_pos = interpolate(e, pos, T.name());
    e_neg = interpolate(e, neg, T.name());

    U_pos = rhoU_pos/rho_pos;
    U_neg = rhoU_neg/rho_neg;

    p_pos = (rho_pos - B0_pos) * rPsi_pos;
    p_neg = (rho_neg - B0_neg) * rPsi_neg;

    phiv_pos = U_pos & mesh.Sf();
    phiv_pos.setOriented(false);
    phiv_neg = U_neg & mesh.Sf();
    phiv_neg.setOriented(false);
}
