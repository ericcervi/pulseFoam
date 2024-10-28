{
    solve
    (
        fvm::ddt(rhoE)
      + fvc::div(phiEp)
    );

    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();

    T = e/Cv;
    T.correctBoundaryConditions();

    e.boundaryFieldRef() = Cv.boundaryField() * T.boundaryField();

    rhoE.boundaryFieldRef() == rho.boundaryField() *
                               (
                                 e.boundaryField() + 0.5*magSqr(U.boundaryField())
                               );
}

