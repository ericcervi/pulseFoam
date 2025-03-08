for (int mI=0; mI<materials; mI++) 
{
    solve(fvm::ddt(rho_m[mI]) + fvc::div(phi));
}

