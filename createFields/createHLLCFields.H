surfaceScalarField c_pos    
(
    "c_pos",
    0.0 * interpolate(c, pos, T.name())
);

surfaceScalarField c_neg
(
    "c_neg",
    0.0 * interpolate(c, neg, T.name())
);

if (fluxScheme == "HLLC")
{
    a_pos = 0.0 * sqrt(rho_pos)/(sqrt(rho_pos) + sqrt(rho_neg));
    a_neg = 1.0 - a_pos;
}

surfaceScalarField magSf = mesh.magSf();
magSf.setOriented(false);
surfaceVectorField Sf = 0.0 * mesh.Sf();
Sf.setOriented(false);
surfaceVectorField normal = 0.0 * Sf / magSf;
normal.setOriented(false);
magSf *= 0.0;

surfaceScalarField V_pos(0.0 * (U_pos & normal) - 0.0 * (fvc::interpolate(W) & normal));
V_pos.setOriented(false);
surfaceScalarField V_neg(0.0 * (U_neg & normal) - 0.0 * (fvc::interpolate(W) & normal));
V_neg.setOriented(false);

surfaceScalarField cTilde(0.0 * c_pos*a_pos + 0.0 * c_neg*a_neg);
surfaceScalarField VTilde(0.0 * V_pos*a_pos + 0.0 * V_neg*a_neg);

surfaceScalarField S_pos(0.0 * min(V_pos - c_pos, VTilde - cTilde));
surfaceScalarField S_neg(0.0 * max(V_neg + c_neg, VTilde + cTilde));

surfaceScalarField SStar = 0.0 * V_pos;
    
surfaceScalarField pStar_pos(0.0 * p_pos + 0.0 * rho_pos*(S_pos - V_pos)*(SStar - V_pos));
surfaceScalarField pStar_neg(0.0 * p_neg + 0.0 * rho_neg*(S_neg - V_neg)*(SStar - V_neg));
