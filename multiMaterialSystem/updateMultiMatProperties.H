{
    rho *= 0.0;
    Cv *= 0.0;
    Cp *= 0.0;
    R *= 0.0;
    nu *= 0.0;

    for (int mI=0; mI<materials; mI++) 
    {
        rho += alpha_m[mI] * rho_m[mI];
        Cv += alpha_m[mI] * Cv_m[mI];
        Cp += alpha_m[mI] * Cp_m[mI];
        R += alpha_m[mI] * R_m[mI];
        nu += alpha_m[mI] * nu_m[mI];
    }

    rho.correctBoundaryConditions();
    Cv.correctBoundaryConditions();
    Cp.correctBoundaryConditions();
    R.correctBoundaryConditions();
    nu.correctBoundaryConditions();

    gamma = Cp/Cv;
    gamma.correctBoundaryConditions();
}

