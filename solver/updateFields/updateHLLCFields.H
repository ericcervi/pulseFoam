{
    c_pos = interpolate(c, pos, T.name());
  
    c_neg = interpolate(c, neg, T.name());
 
    a_pos = sqrt(rho_pos)/(sqrt(rho_pos) + sqrt(rho_neg));
 
    a_neg = 1.0 - a_pos;
}
