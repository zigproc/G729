function [gp_opt,gc_opt] = compute_optimal_codebook_gains(x,x_fcb,y,z)

  % x is target signal (u3 in book)
  % y is zero state response of ACB vector with weighted synthesis filter (WSF) (y2o in book)
  % z is zero state response of FCB vector (with pitch synthesis) with WSF (y1o in book)

  % setup system of equations 
  % M*[gp_opt;gc_opt] = b
  
  gp_opt = 0;
  gc_opt = 0;
  
  M = [ y*y' z*y'; y*z' z*z'];
  b = [ x*y'; x*z'];
  
  if( det(M) ~= 0)
    g = inv(M)*b;
    gp_opt = g(1);
    gc_opt = g(2);
  else  
    gp_opt = 0;
    gc_opt = (x_fcb * z') / (z*z');  % workaround in case no solution to system
  end
  
end