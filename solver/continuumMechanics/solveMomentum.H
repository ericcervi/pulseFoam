{
    solve(fvm::ddt(rhoU) + fvc::div(phiUp) - rho * gravity);

    U.ref() = rhoU()/rho();
    U.correctBoundaryConditions();
    W = U * (1.0 - eulRel * eulerian / max(eulerian));
    W.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

    if (!inviscid)
    {
        solve
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff, U)
          - fvc::div(tauMC)
        );
        rhoU = rho*U;
    }
}
