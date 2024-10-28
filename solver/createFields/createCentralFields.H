surfaceScalarField cSf_pos
(
    "cSf_pos",
    0.0 * interpolate(c, pos, T.name())*mesh.magSf()
);

surfaceScalarField cSf_neg
(
    "cSf_neg",
    0.0 * interpolate(c, neg, T.name())*mesh.magSf()
);

surfaceScalarField ap
(
    "ap",
    0.0 * max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
);

surfaceScalarField am
(
    "am",
    0.0 * min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
);

surfaceScalarField a_pos
(
    IOobject
    (
        "a_pos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

surfaceScalarField amaxSf("amaxSf", 0.0 * max(mag(am), mag(ap)));

surfaceScalarField aSf("aSf", 0.0 * am*a_pos);

if (fluxScheme == "Tadmor")
{
    aSf = -0.5*amaxSf;
    a_pos = 0.5;
}

surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

phiv_pos *= a_pos;
phiv_neg *= a_neg;

surfaceScalarField aphiv_pos("aphiv_pos", 0.0 * phiv_pos);
surfaceScalarField aphiv_neg("aphiv_neg", 0.0 * phiv_neg);

amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
