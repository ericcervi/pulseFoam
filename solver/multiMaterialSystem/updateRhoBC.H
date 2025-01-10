rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField() + B0.boundaryField();

for (int mI=0; mI<materials; mI++) 
{
    if(EOS_m[mI]==0)
    {
        rho_m[mI].boundaryFieldRef() == p.boundaryField() / R_m[mI].value() / T.boundaryField();
    }
    if(EOS_m[mI]==1)
    {
        rho_m[mI].boundaryFieldRef() == (rho.boundaryField()/rho.boundaryField() - cavitated.boundaryField()) *
                                        (
                                          psi_m[mI].value() * p.boundaryField()
                                        - psi_m[mI].value() * p0.value()
                                        - beta_m[mI].value() * (T.boundaryField() - T0.value())
                                        + rho0_m[mI].value()
                                        )
                                        + cavitated.boundaryField() * p.boundaryField() / R_m[mI].value() / T.boundaryField();
    }
}
