forAll(phi, faceI)
{
    if (S_pos[faceI] > 0)
    {
        for (int mI=0; mI<materials; mI++) 
        {
            evaluateLR(alphaPhiv_m[mI][faceI], phi_m[mI][faceI], phiv[faceI], phi[faceI], phiUp[faceI], phiEp[faceI], V_pos[faceI], alpha_pos_m[mI][faceI], rho_pos_m[mI][faceI], rho_pos[faceI], U_pos[faceI], p_pos[faceI], e_pos[faceI], magSf[faceI], Sf[faceI], void_limit, mI);
        }
    }
    else if (SStar[faceI] > 0)
    {
        for (int mI=0; mI<materials; mI++) 
        {
            evaluateStar(alphaPhiv_m[mI][faceI], phi_m[mI][faceI], phiv[faceI], phi[faceI], phiUp[faceI], phiEp[faceI], V_pos[faceI], S_pos[faceI], SStar[faceI], pStar_pos[faceI], normal[faceI], alpha_pos_m[mI][faceI], rho_pos_m[mI][faceI], rho_pos[faceI], U_pos[faceI], p_pos[faceI], e_pos[faceI], magSf[faceI], Sf[faceI], void_limit, mI);    
        }
    }
    else if (S_neg[faceI] > 0)
    {
        for (int mI=0; mI<materials; mI++) 
        {
            evaluateStar(alphaPhiv_m[mI][faceI], phi_m[mI][faceI], phiv[faceI], phi[faceI], phiUp[faceI], phiEp[faceI], V_neg[faceI], S_neg[faceI], SStar[faceI], pStar_neg[faceI], normal[faceI], alpha_neg_m[mI][faceI], rho_neg_m[mI][faceI], rho_neg[faceI], U_neg[faceI], p_neg[faceI], e_neg[faceI], magSf[faceI], Sf[faceI], void_limit, mI);    
        }
    }
    else
    {
        for (int mI=0; mI<materials; mI++) 
        {
            evaluateLR(alphaPhiv_m[mI][faceI], phi_m[mI][faceI], phiv[faceI], phi[faceI], phiUp[faceI], phiEp[faceI], V_neg[faceI], alpha_neg_m[mI][faceI], rho_neg_m[mI][faceI], rho_neg[faceI], U_neg[faceI], p_neg[faceI], e_neg[faceI], magSf[faceI], Sf[faceI], void_limit, mI);
        }
    }
}

forAll(phi.boundaryFieldRef(), patchID)
{
	forAll(phi.boundaryFieldRef()[patchID],faceI) 
    {    
        auto& phiv_pf = phiv.boundaryFieldRef()[patchID][faceI];
        auto& phi_pf = phi.boundaryFieldRef()[patchID][faceI];
        auto& phiUp_pf = phiUp.boundaryFieldRef()[patchID][faceI];
        auto& phiEp_pf = phiEp.boundaryFieldRef()[patchID][faceI];

        if (S_pos.boundaryFieldRef()[patchID][faceI] > 0)
        {
            const auto& rho_pos_pf = rho_pos.boundaryFieldRef()[patchID][faceI];          
            const auto& U_pos_pf = U_pos.boundaryFieldRef()[patchID][faceI]; 
            const auto& V_pos_pf = V_pos.boundaryFieldRef()[patchID][faceI];
            const auto& e_pos_pf = e_pos.boundaryFieldRef()[patchID][faceI];
            const auto& p_pos_pf = p_pos.boundaryFieldRef()[patchID][faceI];
            const auto& Sf_pf = Sf.boundaryFieldRef()[patchID][faceI];
            const auto& magSf_pf = magSf.boundaryFieldRef()[patchID][faceI];

            for (int mI=0; mI<materials; mI++) 
            {
                auto& alpha_pos_pf_m = alpha_pos_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& rho_pos_pf_m = rho_pos_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& alphaPhiv_pf_m = alphaPhiv_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& phi_pf_m = phi_m[mI].boundaryFieldRef()[patchID][faceI];

                evaluateLR(alphaPhiv_pf_m, phi_pf_m, phiv_pf, phi_pf, phiUp_pf, phiEp_pf, V_pos_pf, alpha_pos_pf_m, rho_pos_pf_m, rho_pos_pf, U_pos_pf, p_pos_pf, e_pos_pf, magSf_pf, Sf_pf, void_limit, mI);
            }
        }
        else if (SStar.boundaryFieldRef()[patchID][faceI] > 0)
        {            
            const auto& rho_pos_pf = rho_pos.boundaryFieldRef()[patchID][faceI];        
            const auto& U_pos_pf = U_pos.boundaryFieldRef()[patchID][faceI]; 
            const auto& V_pos_pf = V_pos.boundaryFieldRef()[patchID][faceI];
            const auto& e_pos_pf = e_pos.boundaryFieldRef()[patchID][faceI];
            const auto& p_pos_pf = p_pos.boundaryFieldRef()[patchID][faceI];
            const auto& Sf_pf = Sf.boundaryFieldRef()[patchID][faceI];
            const auto& magSf_pf = magSf.boundaryFieldRef()[patchID][faceI];        
            const auto& normal_pf = normal.boundaryFieldRef()[patchID][faceI];

            const auto& S_pos_pf = S_pos.boundaryFieldRef()[patchID][faceI];        
            const auto& SStar_pf = SStar.boundaryFieldRef()[patchID][faceI];        
            const auto& pStar_pos_pf = pStar_pos.boundaryFieldRef()[patchID][faceI];        
                        
            for (int mI=0; mI<materials; mI++) 
            {
                auto& alpha_pos_pf_m = alpha_pos_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& rho_pos_pf_m = rho_pos_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& alphaPhiv_pf_m = alphaPhiv_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& phi_pf_m = phi_m[mI].boundaryFieldRef()[patchID][faceI];
            
                evaluateStar(alphaPhiv_pf_m, phi_pf_m, phiv_pf, phi_pf, phiUp_pf, phiEp_pf, V_pos_pf, S_pos_pf, SStar_pf, pStar_pos_pf, normal_pf, alpha_pos_pf_m, rho_pos_pf_m, rho_pos_pf, U_pos_pf, p_pos_pf, e_pos_pf, magSf_pf, Sf_pf, void_limit, mI);
            }
        }
        else if (S_neg.boundaryFieldRef()[patchID][faceI] > 0)
        {
            const auto& rho_neg_pf = rho_neg.boundaryFieldRef()[patchID][faceI];        
            const auto& U_neg_pf = U_neg.boundaryFieldRef()[patchID][faceI]; 
            const auto& V_neg_pf = V_neg.boundaryFieldRef()[patchID][faceI];
            const auto& e_neg_pf = e_neg.boundaryFieldRef()[patchID][faceI];
            const auto& p_neg_pf = p_neg.boundaryFieldRef()[patchID][faceI];
            const auto& Sf_pf = Sf.boundaryFieldRef()[patchID][faceI];
            const auto& magSf_pf = magSf.boundaryFieldRef()[patchID][faceI];        
            const auto& normal_pf = normal.boundaryFieldRef()[patchID][faceI];        
            
            const auto& S_neg_pf = S_neg.boundaryFieldRef()[patchID][faceI];        
            const auto& SStar_pf = SStar.boundaryFieldRef()[patchID][faceI];        
            const auto& pStar_neg_pf = pStar_neg.boundaryFieldRef()[patchID][faceI];        

            for (int mI=0; mI<materials; mI++) 
            {
                auto& alpha_neg_pf_m = alpha_neg_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& rho_neg_pf_m = rho_neg_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& alphaPhiv_pf_m = alphaPhiv_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& phi_pf_m = phi_m[mI].boundaryFieldRef()[patchID][faceI];
            
                evaluateStar(alphaPhiv_pf_m, phi_pf_m, phiv_pf, phi_pf, phiUp_pf, phiEp_pf, V_neg_pf, S_neg_pf, SStar_pf, pStar_neg_pf, normal_pf, alpha_neg_pf_m, rho_neg_pf_m, rho_neg_pf, U_neg_pf, p_neg_pf, e_neg_pf, magSf_pf, Sf_pf, void_limit, mI);
            }
        }
        else
        {      
            const auto& rho_neg_pf = rho_neg.boundaryFieldRef()[patchID][faceI];          
            const auto& U_neg_pf = U_neg.boundaryFieldRef()[patchID][faceI]; 
            const auto& V_neg_pf = V_neg.boundaryFieldRef()[patchID][faceI];
            const auto& e_neg_pf = e_neg.boundaryFieldRef()[patchID][faceI];
            const auto& p_neg_pf = p_neg.boundaryFieldRef()[patchID][faceI];
            const auto& Sf_pf = Sf.boundaryFieldRef()[patchID][faceI];
            const auto& magSf_pf = magSf.boundaryFieldRef()[patchID][faceI];

            for (int mI=0; mI<materials; mI++) 
            {
                auto& alpha_neg_pf_m = alpha_neg_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& rho_neg_pf_m = rho_neg_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& alphaPhiv_pf_m = alphaPhiv_m[mI].boundaryFieldRef()[patchID][faceI];
                auto& phi_pf_m = phi_m[mI].boundaryFieldRef()[patchID][faceI];
            
                evaluateLR(alphaPhiv_pf_m, phi_pf_m, phiv_pf, phi_pf, phiUp_pf, phiEp_pf, V_neg_pf, alpha_neg_pf_m, rho_neg_pf_m, rho_neg_pf, U_neg_pf, p_neg_pf, e_neg_pf, magSf_pf, Sf_pf, void_limit, mI);
            }
        }
    }
}

// Make relative to mesh motion
if (eulRel<1.0) 
{  
    surfaceScalarField meshPhi = fvc::interpolate(W) & mesh.Sf();
    meshPhi.setOriented(false);       
    phiv += meshPhi;
}
