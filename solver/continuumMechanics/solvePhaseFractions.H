{
    if (materials > 1)
    {
        for (int mI=0; mI<materials; mI++) 
        {
            solve(fvm::ddt(alpha_m[mI]) + fvc::div(alphaPhiv_m[mI]) - alpha_m[mI] * fvc::div(phiv));
        }
    }
}

