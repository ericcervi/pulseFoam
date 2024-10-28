{
    magSf = mesh.magSf();
    magSf.setOriented(false);
    Sf = mesh.Sf();
    Sf.setOriented(false);
    normal = Sf/magSf;
    normal.setOriented(false);

    V_pos = (U_pos & normal) - (fvc::interpolate(W) & normal);
    V_pos.setOriented(false);
    V_neg = (U_neg & normal) - (fvc::interpolate(W) & normal);
    V_neg.setOriented(false);
       
    cTilde = c_pos*a_pos + c_neg*a_neg;
    VTilde = V_pos*a_pos + V_neg*a_neg;
        
    S_pos = min(V_pos - c_pos, VTilde - cTilde);
    S_neg = max(V_neg + c_neg, VTilde + cTilde);
    
    SStar =
    (
        p_neg - p_pos
      + rho_pos*V_pos*(S_pos - V_pos)
      - rho_neg*V_neg*(S_neg - V_neg)
    )
   /(rho_pos*(S_pos - V_pos) - rho_neg*(S_neg - V_neg));

    
    pStar_pos = p_pos + rho_pos*(S_pos - V_pos)*(SStar - V_pos);
    pStar_neg = p_neg + rho_neg*(S_neg - V_neg)*(SStar - V_neg);
}
