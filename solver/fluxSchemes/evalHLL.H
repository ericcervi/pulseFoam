namespace Foam
{
    void evaluateCenter(scalar& alphaPhiv_m, scalar& phi_m, scalar& phiv, scalar& phi, vector& phiUp, scalar& phiEp, scalar V_neg, scalar V_pos, scalar S_neg, scalar S_pos, vector normal, scalar alpha_neg_m, scalar alpha_pos_m, scalar rho_neg_m, scalar rho_pos_m, scalar rho_neg, scalar rho_pos, vector U_neg, vector U_pos, scalar p_neg, scalar p_pos, scalar e_neg, scalar e_pos, scalar magSf, vector Sf, scalar void_limit, int mI)
    {
        scalar void_corr_pos = 1.0;
        scalar void_corr_neg = 1.0;
        if (rho_pos < void_limit)
        {
            void_corr_pos = 0.0;
        }
        if (rho_neg < void_limit)
        {
            void_corr_neg = 0.0;
        }

        alphaPhiv_m = (  
                  S_neg * alpha_pos_m * (V_pos * magSf) * void_corr_neg
                - S_pos * alpha_neg_m * (V_neg * magSf) * void_corr_pos
                + S_pos * S_neg * (alpha_neg_m * void_corr_neg - alpha_pos_m * void_corr_pos) * magSf
                ) / (S_neg - S_pos);
                
        phi_m = (  
                  S_neg * rho_pos_m * (V_pos * magSf) * void_corr_neg
                - S_pos * rho_neg_m * (V_neg * magSf) * void_corr_pos
                + S_pos * S_neg * (rho_neg_m * void_corr_neg - rho_pos_m * void_corr_pos) * magSf
                ) / (S_neg - S_pos);

        if (mI==0)
        {
            phiv = (  
                     S_neg * (V_pos * magSf) * void_corr_neg
                    - S_pos * (V_neg * magSf) * void_corr_pos
                   ) / (S_neg - S_pos);

            phi = (  
                    S_neg * rho_pos * (V_pos * magSf) * void_corr_neg
                  - S_pos * rho_neg * (V_neg * magSf) * void_corr_pos
                  + S_pos * S_neg * (rho_neg * void_corr_neg - rho_pos * void_corr_pos) * magSf
                  ) / (S_neg - S_pos);

            phiUp = (  
                      S_neg * (rho_pos * U_pos * (V_pos * magSf) + p_pos * Sf) * void_corr_neg
                    - S_pos * (rho_neg * U_neg * (V_neg * magSf) + p_neg * Sf) * void_corr_pos
                    + S_pos * S_neg * (rho_neg * U_neg * void_corr_neg - rho_pos * U_pos * void_corr_pos) * magSf
                    ) / (S_neg - S_pos);
            
            phiEp = (  
                      S_neg * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) * (V_pos * magSf) * void_corr_neg
                    - S_pos * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) * (V_neg * magSf) * void_corr_pos
                    + S_pos * S_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) * void_corr_neg - rho_pos * (e_pos + 0.5 * magSqr(U_pos)) * void_corr_pos) * magSf
                    ) / (S_neg - S_pos);
        }  
    }  
}
