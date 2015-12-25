function T = open_loop_pitch_analysis(sw,ZI)
  
  assert(length(sw) == (80+143));
  
  R1 = zeros(1,length(80:143));
  R2 = zeros(1,length(40:79));
  R3 = zeros(1,length(20:39));
  
  % [sw(0)...sw(79)]
  sw_fixed = sw(ZI:ZI+80-1);
  
  % Equation (34) 
  % compute for ranges [80,143], [40,79], [20,39]
  for k=80:143
    R1(k-80+1) = sw_fixed * sw(ZI-k:ZI-k+80-1)';
  end
  
  for k=40:79
    R2(k-40+1) = sw_fixed * sw(ZI-k:ZI-k+80-1)';
  end 
 
  for k=20:39
    R3(k-20+1) = sw_fixed * sw(ZI-k:ZI-k+80-1)';
  end  
  
  [R1_max,t1] = max(R1);
  [R2_max,t2] = max(R2);
  [R3_max,t3] = max(R3);
  
  % convert to actual values, based on starting value of range
  t1 = t1 - 1 + 80;
  t2 = t2 - 1 + 40;
  t3 = t3 - 1 + 20;
  
  % normalize
  % Equation (35)
  sw_t1  = sw(ZI-t1:ZI-t1+80-1);
  D      = sw_t1 * sw_t1';
  if( D ~= 0)
    R1_max = R1_max / sqrt(D);
  end
  
  sw_t2  = sw(ZI-t2:ZI-t2+80-1);
  D      = sw_t2 * sw_t2';
  if( D ~= 0)
    R2_max = R2_max / sqrt(D);  
  end
  
  sw_t3  = sw(ZI-t3:ZI-t3+80-1);
  D      = sw_t3 * sw_t3';
  if( D ~= 0)
    R3_max = R3_max / sqrt(D); 
  end
  % compute winner, using weighted values to favor lower delay
  % values
  T = t1;
  R_max = R1_max;

  if( R2_max >= 0.85*R_max)
    R_max = R2_max;
    T = t2;
  end

  if( R3_max >= 0.85*R_max)
    R_max = R3_max;
    T = t3;
  end  
    
  
end