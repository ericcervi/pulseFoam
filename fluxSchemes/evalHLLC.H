namespace Foam
{
    void evaluateLR(scalar& alphaPhiv_m, scalar& phi_m, scalar& phiv, scalar& phi, vector& phiUp, scalar& phiEp, scalar V_dir, scalar alpha_dir_m, scalar rho_dir_m, scalar rho_dir, vector U_dir, scalar p_dir, scalar e_dir, scalar magSf, vector Sf, scalar void_limit, int mI)
    {
        scalar void_corr = 1.0;
        if (rho_dir < void_limit)
        {
            void_corr = 0.0;
        }

        alphaPhiv_m = V_dir * magSf * alpha_dir_m;
        phi_m = rho_dir_m * V_dir * magSf * void_corr;
            
        if (mI==0)
        {       
            phiv = V_dir * magSf;
            phi = rho_dir * phiv * void_corr;

            phiUp = (phi * U_dir + p_dir * Sf) * void_corr;
            phiEp = (phi * (e_dir + 0.5 * magSqr(U_dir)) + p_dir * V_dir * magSf) * void_corr;
        }
    }
    
    void evaluateStar(scalar& alphaPhiv_m, scalar& phi_m, scalar& phiv, scalar& phi, vector& phiUp, scalar& phiEp, scalar V_dir, scalar S_dir, scalar SStar, scalar pStar_dir, vector normal, scalar alpha_dir_m, scalar rho_dir_m, scalar rho_dir, vector U_dir, scalar p_dir, scalar e_dir, scalar magSf, vector Sf, scalar void_limit, int mI)
    {
        scalar void_corr = 1.0;
        if (rho_dir < void_limit)
        {
            void_corr = 0.0;
        }

        scalar f = (S_dir - V_dir) / (S_dir - SStar);
        alphaPhiv_m = SStar * f * magSf * alpha_dir_m;
        phi_m = rho_dir_m * SStar * f * magSf * void_corr;

        if (mI==0)
        {
            phiv = SStar * f * magSf;
            phi = rho_dir * phiv * void_corr;

            vector rhoUStar = f * rho_dir * (U_dir - (V_dir - SStar) * normal) * magSf;

            scalar rhoEStar = f * rho_dir * (e_dir + 0.5 * magSqr(U_dir)) * magSf
                            + (pStar_dir * SStar - p_dir * V_dir) / (S_dir - V_dir) * magSf;

            phiUp = rho_dir * U_dir * (V_dir * magSf)
                  + p_dir * Sf
                  + S_dir * (rhoUStar - rho_dir * U_dir * magSf);

            phiUp *= void_corr;

            phiEp = rho_dir * (e_dir + 0.5 * magSqr(U_dir)) * (V_dir * magSf)
                  + p_dir * (V_dir * magSf)
                  + S_dir * (rhoEStar - rho_dir * (e_dir + 0.5 * magSqr(U_dir)) * magSf);

            phiEp *= void_corr;
        }  
    }  
}
