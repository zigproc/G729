function v = compute_adaptive_codebook_vector(u,ZI,T_int,T_frac,interp_b)

  % 
  assert(length(u) == (143 + 40 + 9));
  assert(length(interp_b) == 60);
  assert(any(T_frac == -1:1));
  
  % polyphase components of interpolation filter for R(k)_t
  % Equation (39)
  b0 = fliplr(interp_b(1:3:end));  % t = 0 
  b1 = fliplr(interp_b(2:3:end));  % t = 1 (frac 1/3)
  b2 = fliplr(interp_b(3:3:end));  % t = 2 (frac 2/3)
  
  % convert fractional pitch delay, encoded as (T_int,T_frac), where T_frac is -1,0,1
  % To compute adaptive codebook vector we need to access u(n-k) where k is
  % fractional resolution.  So we need an integer number of samples in the past, 
  % moved up by a positive fractional amount 0,1/3,or 2/3
  
  % Encoding -> Value represented -> how to index as u(n-k+t) where k int, t frac
  % 59,-1 -> 58 2/3 -> -59 + 1/3 :    u(n-58.66666) = u(n-59+0.3333)
  % 59,0 -> 59 -> -59
  % 59,1 -> 59 1/3 -> -60 + 2/3
  
  if( T_frac == -1)
    T_frac = 1;
  elseif( T_frac == 1)
    T_frac = 2;
    T_int = T_int + 1;
  end
  
  % select which filter we need to use for interpolation
  switch(T_frac)
    case 0
      b = b0;
    case 1
      b = b1;
    case 2
      b = b2;
    otherwise 
      assert(false);
  endswitch
  
  % clear out the "new" samples
  % we will interpolate in-place, since for small pitch delays we will need
  % to use newly interpolated samples in this range
  u(ZI:ZI+40-1) = zeros(1,40);
  
  % compute "new" adaptive codebook vector samples by interpolating past
  % excitation
  
  assert(ZI-T_int-9 > 0)
  
  for n=0:39
    % u(n-k-9)...u(n-k+10), where k is integer valued lag
    s = u(ZI+n-T_int-9:ZI+n-T_int+10);  % 20 samples surrounding integer lag
    u(ZI+n) = s*b';                     % interpolated sample value
  end
  
  % the adaptive codebook vector is the new samples
  v = u(ZI:ZI+40-1);
end