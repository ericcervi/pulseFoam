{
    solve
    (
        fvm::ddt(rhoE)
      + fvc::div(phiEp)
    );

    e = rhoE/rho - 0.5*magSqr(U);

    /*forAll(e, cellI)
    {
        if (e[cellI] <= 0.0)
        {
            e[cellI] = 1e-9;
        }
    }*/
    e.correctBoundaryConditions();

    T = e/Cv;
    T.correctBoundaryConditions();

    e.boundaryFieldRef() = Cv.boundaryField() * T.boundaryField();

    rhoE.boundaryFieldRef() == rho.boundaryField() *
                               (
                                 e.boundaryField() + 0.5*magSqr(U.boundaryField())
                               );
}

