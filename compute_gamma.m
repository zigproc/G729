function [gamma1,gamma2,flat] = compute_gamma(o,prev_flat,w)

  assert(length(o) == 2);
  assert(length(w) == 10);
    
  % Equation (30)
  if( prev_flat == 1 && o(1) < -1.74 && o(2) > 0.65 )
      flat = 0;
  elseif( prev_flat == 0 && (o(1) > -1.52 || o(2) < 0.43))
      flat = 1;
  else
      flat = prev_flat;
  end
  
  if( flat == 1)
    gamma1 = 0.94;
    gamma2 = 0.6;
  else
    gamma1 = 0.98;
    % Equation (31)
    dmin = min(w(2:10)-w(1:9));
    % Equation (32)
    gamma2 = -6.0*dmin + 1.0;
    gamma2 = max(gamma2,0.4);
    gamma2 = min(gamma2,0.7);
  end
    
end
