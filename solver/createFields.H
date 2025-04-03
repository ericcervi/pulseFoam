#include "createRDeltaT.H"

int materials=runTime.controlDict().lookupOrDefault<scalar>("materials", 1);

IOdictionary materialProperties
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Foam::dimensionedVector gravity = materialProperties.lookupOrDefault<Foam::dimensionedVector>(
    "gravity", Foam::dimensionedVector("gravity", Foam::dimensionSet(0, 1, -2, 0, 0, 0, 0), Foam::vector(0.0, 0.0, 0.0)));

std::vector<int> EOS_m{};
std::vector<dimensionedScalar> R_m{};
std::vector<dimensionedScalar> beta_m{};
std::vector<dimensionedScalar> psi_m{};
std::vector<dimensionedScalar> Cv_m{};
std::vector<dimensionedScalar> Cp_m{};
std::vector<dimensionedScalar> rho0_m{};
std::vector<dimensionedScalar> nu_m{};

// multi-material simulation control options
scalar void_limit = materialProperties.lookupOrDefault<scalar>("void_limit", -1.0);
scalar p_cavitation = materialProperties.lookupOrDefault<scalar>("p_cavitation", -1.0e30);
int cavitating_id = materialProperties.lookupOrDefault<int>("cavitating_id", 1000);
bool boiling = materialProperties.lookupOrDefault<bool>("boiling", false);

// multi-material properties
for (int mI=0; mI<materials; mI++) 
{
    scalar EOS_value = materialProperties.lookupOrDefault<scalar>("EOS_m" + Foam::name(mI), 0.0);
    EOS_m.push_back(EOS_value);

    dimensionedScalar R_value(materialProperties.lookup("R_m" + Foam::name(mI)));
    R_m.push_back(R_value);

    dimensionedScalar nu_value(materialProperties.lookupOrDefault<dimensionedScalar>(
        "nu_m" + Foam::name(mI), dimensionedScalar("nu_value", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ));
    nu_m.push_back(nu_value);

    if (EOS_value == 1)
    {
	    dimensionedScalar beta_value(materialProperties.lookup("beta_m" + Foam::name(mI)));
        beta_m.push_back(beta_value);
    
	    dimensionedScalar psi_value(materialProperties.lookup("psi_m" + Foam::name(mI)));
        psi_m.push_back(psi_value);

	    dimensionedScalar rho0_value(materialProperties.lookup("rho0_m" + Foam::name(mI)));
        rho0_m.push_back(rho0_value);
    }
    else
    { 
        beta_m.push_back(0.0);
        psi_m.push_back(0.0);
        rho0_m.push_back(0.0);
    }

	dimensionedScalar Cv_value(materialProperties.lookup("Cv_m" + Foam::name(mI)));
    Cv_m.push_back(Cv_value);
    
	dimensionedScalar Cp_value(materialProperties.lookup("Cp_m" + Foam::name(mI)));
    Cp_m.push_back(Cp_value);
}

dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar T0(materialProperties.lookup("T0"));

PtrList<volScalarField> alpha_m(materials);
forAll(alpha_m, mI)
{
	alpha_m.set
	(
		mI,
		new volScalarField
		(
	        IOobject
            (
	            "alpha_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh
		)
	);
}

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Acceleration field
// Implemented to address ANS reviewer comment
// To be removed
/*volVectorField U_old = U;

volVectorField a
(
    IOobject
    (
        "a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (U - U_old) / runTime.deltaT()
);*/

PtrList<volScalarField> rho_m(materials);
forAll(rho_m, mI)
{
	rho_m.set
	(
		mI,
		new volScalarField
		(
	        IOobject
            (
	            "rho_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh
		)
	);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha_m[0]*rho_m[0]*0.0
);

volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha_m[0]*nu_m[0]
);

nu *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    nu += alpha_m[0]*nu_m[0];
}
nu.correctBoundaryConditions();

bool inviscid = true;
if (max(nu).value() > 0.0)
{
    inviscid = false;
}

volScalarField powDens
(
    IOobject
    (
        "powDens",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("powDens",dimensionSet(1,-1,-3,0,0,0,0),0)
);

volScalarField cavitated
(
    IOobject
    (
        "cavitated",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha_m[0]*rho_m[0]*0.0
);

volScalarField boiled
(
    IOobject
    (
        "boiled",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha_m[0]*rho_m[0]*0.0
);

rho *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    rho += alpha_m[0]*rho_m[0];
}
rho.correctBoundaryConditions();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha_m[0]*R_m[0]*0.0
);

R *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    R += alpha_m[0]*R_m[0];
}
R.correctBoundaryConditions();

volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha_m[0]*Cv_m[0]*0.0
);

Cv *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    Cv += alpha_m[0]*Cv_m[0];
}
Cv.correctBoundaryConditions();

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha_m[0]*Cp_m[0]*0.0
);

Cp *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    Cp += alpha_m[0]*Cp_m[0];
}
Cp.correctBoundaryConditions();

//// to be removed
rho *= 0.0;
Cv *= 0.0;
Cp *= 0.0;
R *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    rho += alpha_m[mI] * rho_m[mI];
    Cv += alpha_m[mI] * Cv_m[mI];
    Cp += alpha_m[mI] * Cp_m[mI];
    R += alpha_m[mI] * R_m[mI];
}
rho.correctBoundaryConditions();
Cv.correctBoundaryConditions();
Cp.correctBoundaryConditions();
R.correctBoundaryConditions();
////

volScalarField gamma = Cp/Cv;
gamma.correctBoundaryConditions();

volScalarField e = Cv*T + 0.0*p/rho;
e.correctBoundaryConditions();

volScalarField B0 = rho * 0.0;
B0.correctBoundaryConditions();

volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("psi",dimensionSet(0,-2,2,0,0,0,0),0)
);

psi *= 0.0;
B0 *= 0.0;
for (int mI=0; mI<materials; mI++) 
{
    if(EOS_m[mI]==0)
    {
        psi += alpha_m[mI] / R_m[mI] / T;
    }
    if(EOS_m[mI]==1)
    {
        psi += alpha_m[mI] * psi_m[mI];
        B0 +=  alpha_m[mI] * (rho0_m[mI] - beta_m[mI]*(T - T0) - psi_m[mI]*p0);
    }
}

forAll(psi, cellI)
{
    if (psi[cellI] <= 0.0)
    {
        psi[cellI] = 1e-9;
    }
}
psi.correctBoundaryConditions();
B0.correctBoundaryConditions();

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

PtrList<surfaceScalarField> alpha_pos_m(materials);
forAll(alpha_pos_m, mI)
{
	alpha_pos_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "alpha_pos_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            interpolate(alpha_m[0], pos)
		)
	);
}

PtrList<surfaceScalarField> alpha_neg_m(materials);
forAll(alpha_neg_m, mI)
{
	alpha_neg_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "alpha_neg_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            interpolate(alpha_m[0], neg)
		)
	);
}

PtrList<surfaceScalarField> rho_pos_m(materials);
forAll(rho_pos_m, mI)
{
	rho_pos_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "rho_pos_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            interpolate(rho_m[0], pos)
		)
	);
}

PtrList<surfaceScalarField> rho_neg_m(materials);
forAll(rho_neg_m, mI)
{
	rho_neg_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "rho_neg_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            interpolate(rho_m[0], neg)
		)
	);
}

PtrList<surfaceScalarField> alphaPhiv_m(materials);
forAll(alphaPhiv_m, mI)
{
	alphaPhiv_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "alphaPhiv_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            fvc::flux(U)*0.0
		)
	);
}

PtrList<surfaceScalarField> phi_m(materials);
forAll(phi_m, mI)
{
	phi_m.set
	(
		mI,
		new surfaceScalarField
		(
	        IOobject
            (
	            "phi_m"+Foam::name(mI),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            fvc::flux(rhoU)*0.0
		)
	);
}

surfaceScalarField phi("phi", fvc::flux(rhoU));

scalar eulRel = runTime.controlDict().lookupOrDefault<scalar>("eulRel", 0);
scalar externalCycles = runTime.controlDict().lookupOrDefault<scalar>("externalCycles", 1);

volScalarField hllc_frac
(
    IOobject
    (
        "hllc_frac",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hllc_frac", dimless, 1.0)
);


volScalarField eulerian
(
    IOobject
    (
        "eulerian",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField W = U * (1.0 - eulRel * eulerian / max(eulerian));
W.correctBoundaryConditions();

pointMesh pMesh(mesh);

pointVectorField pointDisp
(
IOobject
(
        "pointDisp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
),
pMesh,
dimensionSet(0,1,0,0,0,0,0)
); 

pointVectorField pointU
(
IOobject
(
        "pointU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
),
pMesh,
dimensionSet(0,1,-1,0,0,0,0)
); 

pointScalarField pointE
(
IOobject
(
        "pointE",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
),
pMesh,
dimensionSet(0,-2,-1,0,0,0,0)
);

volPointInterpolation pointInterpolation(mesh);

vectorField initialPoints = mesh.points();

vectorField newPoints = mesh.points();

