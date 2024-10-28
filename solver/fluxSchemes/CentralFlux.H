{
    for (int mI=0; mI<materials; mI++) 
    {
        if (materials > 1)
        {
            alphaPhiv_m[mI] = aphiv_pos*alpha_pos_m[mI] + aphiv_neg*alpha_neg_m[mI];
        }
        phi_m[mI] = aphiv_pos*rho_pos_m[mI] + aphiv_neg*rho_neg_m[mI];
    }

    phiv = aphiv_pos + aphiv_neg;

    if (eulRel<1.0) 
    {  
        surfaceScalarField meshPhi = fvc::interpolate(W) & mesh.Sf();
        meshPhi.setOriented(false);       
        phiv += meshPhi;
    }

    phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

    phiU = aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg;
    phiU.setOriented(true);

    phiUp = phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf();

    phiEp =
    (
        "phiEp",
        aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
      + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
      + aSf*p_pos - aSf*p_neg
    );
}
