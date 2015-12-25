function MSE = gain_quantization_calc_MSE(x,y,z,gp,gc)

   % need to minimize Equation (63)
   Ex  = x*x';   
   Ey  = y*y';
   Ez  = z*z';
   Exy = x*y';
   Exz = x*z';
   Eyz = y*z';
   
   MSE = Ex +...
         (gp^2)*Ey  +...
         (gc^2)*Ez  -...
         (2*gp)*Exy -...
         (2*gc)*Exz +...
         (2*gp*gc)*Eyz;

end