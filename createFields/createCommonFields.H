surfaceScalarField rho_pos(interpolate(rho, pos));
surfaceScalarField rho_neg(interpolate(rho, neg));
        
surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

volScalarField rPsi("rPsi", 1.0/psi);
surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

surfaceScalarField B0_pos(interpolate(B0, pos, rho.name()));
surfaceScalarField B0_neg(interpolate(B0, neg, rho.name()));

surfaceScalarField e_pos(interpolate(e, pos, T.name()));
surfaceScalarField e_neg(interpolate(e, neg, T.name()));

surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

surfaceScalarField p_pos("p_pos", (rho_pos - B0_pos) * rPsi_pos);
surfaceScalarField p_neg("p_neg", (rho_neg - B0_neg) * rPsi_neg);

surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
phiv_pos.setOriented(false);
surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
phiv_neg.setOriented(false);

volScalarField c("c", sqrt(gamma*rPsi));
