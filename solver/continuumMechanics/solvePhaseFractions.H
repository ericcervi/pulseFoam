{
    if (materials > 1)
    {
        for (int mI=0; mI<materials; mI++) 
        {
            solve(fvm::ddt(alpha_m[mI]) + fvc::div(alphaPhiv_m[mI]) - alpha_m[mI] * fvc::div(phiv));

            forAll(alpha_m[mI], cellI)
            {
                if (alpha_m[mI][cellI] < 0.0)
                {
                 //   alpha_m[mI][cellI] = 1.0e-6;
                }
            }
            alpha_m[mI].correctBoundaryConditions();
        }
    }
}

