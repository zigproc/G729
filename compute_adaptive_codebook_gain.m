function [gp,y] = compute_adaptive_codebook_gain(v,x,h)

  assert(length(v) == 40);
  assert(length(x) == 40);
  assert(length(h) == 40);

  % convolve v(n) with h(n)
  % Equation (44)  
  y = conv_zero_state(v,h);

  gp = 0;
  D  = y*y';
  % Equation (43)
  
  if(D ~= 0)
    gp = x*y' / (y*y');
  end
  gp = max(0,gp);
  gp = min(gp,1.2);
 
end