{
    solve(fvm::ddt(rhoU) + fvc::div(phiUp) - rho * gravity);

    U.ref() = rhoU()/rho();
    U.correctBoundaryConditions();
    W = U * (1.0 - eulRel * eulerian / max(eulerian));
    W.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
}
