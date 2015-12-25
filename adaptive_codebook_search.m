function [T_int,T_frac] = adaptive_codebook_search(T,subframe_idx,x,u,ZI_U,h,interp_r)

  % TBD.  
  assert(length(u) == (143 + 40 + 9));
  assert(length(interp_r) == 24);
  
  % polyphase components of interpolation filter for R(k)_t
  % Equation (39)
  b0 = fliplr(interp_r(1:3:end));  % t = 0 
  b1 = fliplr(interp_r(2:3:end));  % t = 1 (frac 1/3)
  b2 = fliplr(interp_r(3:3:end));  % t = 2 (frac 2/3)
  
  % set closed loop search range
  if( subframe_idx == 1)
    LL = 3;
    UL = 6;
  else
    LL = 5;
    UL = 9;
  end


  tmin = T - LL;
  if( tmin < 20)
    tmin = 20;
  end
    
  tmax = tmin + UL;
  if( tmax > 143)
    tmax = 143;
    tmin = tmax - UL;
  end
  
  % range to compute R(k), need to include extra samples for the interpolation
  k_min = tmin - 4;
  k_max = tmax + 4;
  k_max = min(k_max,143);


  % R(1) is for tmin-4
  % R(end) is for tmax + 4  
  R = zeros(1,(k_max - k_min + 1));
  yk = zeros(1,40);
 
  % Equation (38)
  % calculate for minimum value
  %uk = u(ZI_U - k_min - 40 + 1  : ZI_U - k_min);
  uk = u(ZI_U - k_min : ZI_U - k_min + 40 - 1);  
  yk = conv_zero_state(h,uk);
  
  R(0+(1)) = compute_R(x,yk);
  
  % recursive relation 
  % y_k(n) = y_{k-1}(n-1) + u(-k)*h(n)
  for k=(k_min+1):k_max
    yk = [0 yk(1:(end-1))] + u(ZI_U - k)*h;
    R(k-k_min+(1)) = compute_R(x,yk);
  end
    
  % find maximum R value over range tmin:tmax.  
  k_idx = (tmin:tmax);;
  [m,idx] = max(R(k_idx - k_min + (1)));
  
  T_int = k_idx(idx);    % integer pitch delay corresponding to max of R in [tmin,tmax]
  T_frac = 0;
  
  % if optimum is < 85, do fractional search
  if( T_int < 85)
    R_frac = zeros(1,5);   % fractional range (-2:2)/3
    k_idx  = T_int-4:T_int+4;         % interpolation for R(k) needs R(k-3) thru R(k+4).
                                      % need max range for k=T_int-1 and k=Tint
    R_vals = R(k_idx - k_min + (1));  % R values needed for interpolation
    
    s0 = R_vals(1:end-1);   % range needed for R(T_int-1) interpolation + 1/3, 2/3
    s1 = R_vals(2:end);     % range needed for R(T_Int) interpolation   + 0, 1/3, 2/3
    
    % Equation (39)
    R_frac(1) = s0 * b1';   % -2/3 = decrement Tint by 1, and add 1/3 frac
    R_frac(2) = s0 * b2';   % -1/3 = decrement Tint by 1, and add 2/3 frac
    R_frac(3) = s1 * b0';   % 0
    R_frac(4) = s1 * b1';   % 1/3
    R_frac(5) = s1 * b2';   % 2/3

    % fractional representation of each index in R_frac array 
    % (1/3 units) 
    % e.g. R_frac(1) corresponds to fraction of -2/3   
    R_frac_vals = -2:2;
    
    [m,idx] = max(R_frac);
    
    T_frac = R_frac_vals(idx);
    
    % adjust so T_frac is -1,0,1
    if( T_frac == -2)
      T_int  -= 1;
      T_frac += 3;
    elseif( T_frac == 2)
      T_int  += 1;
      T_frac -= 3;
    end
   
  end 
end

% Equation (37)
function y = compute_R(x,yk)
  D = sqrt(yk*yk');
  if( D == 0)
    y = 0;
  else
    y = (x * yk') / D;
  end
end
