{
    cSf_pos = interpolate(c, pos, T.name())*mesh.magSf();

    cSf_neg = interpolate(c, neg, T.name())*mesh.magSf();

    ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);

    am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

    a_pos = ap/(ap - am);

    amaxSf = max(mag(am), mag(ap));

    aSf = am*a_pos;

    if (fluxScheme == "Tadmor")
    {
        aSf = -0.5*amaxSf;
        a_pos = 0.5;
    }

    a_neg = 1.0 - a_pos;

    phiv_pos *= a_pos;
    phiv_neg *= a_neg;

    aphiv_pos = phiv_pos - aSf;
    aphiv_neg = phiv_neg + aSf;

    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
}
