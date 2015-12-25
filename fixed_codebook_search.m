function [iPulse,sPulse,c,L4_cnt] = fixed_codebook_search(h,T,x,gp_prev,L4_max)

  assert(length(h) == 40);
  assert(length(x) == 40);
  
  % output contains positions and signs of 4 pulses
  iPulse = zeros(1,4);
  sPulse = zeros(1,4);
  
  B = 0;
  
  if( T < 40)
    B = max(gp_prev,0.2);
    B = min(B,0.8);
    
    N = 40 - T;
    
    % Equation (49)
    % add delayed version of h to last N samples
    h(end-N+1:end) = h(end-N+1:end) + B*h(0+(1):N-1+(1));
  end
  
  % Equation (51)
  % compute convolution matrix H
  % compute Phi(i,j) = H'*H
  hMat = tril(toeplitz(h));
  phi  = hMat'*hMat;
  
  % compute correlation 
  % Equation (52)
  d = compute_correlation(x,h);
  
  dMag  = abs(d);
  
  dSign = ones(1,length(d));
  dSign( d < 0) = -1; 
  
  % Modify the matrix Phi with sign information
  % Equation (56)
  % TBD: spec and book disagree.  Either scale this by 0.5, or scale diag entries
  for i=0:39
    for j=i+1:39
      phi(i+(1),j+(1)) = dSign(i+(1))*dSign(j+(1))*phi(i+(1),j+(1));
    end
  end
  
  % keep matrix symmetric.  need same sign modifications as above
  for i=0:39
    for j=i-1:-1:0
      phi(i+(1),j+(1)) = phi(j+(1),i+(1));   % M(i,j) = M(j,i)
    end
  end  
  
  % Equation (57), scale diagonal entries by 2
  % TBD: spec and book disagree.  Scale this by 0.5 or scale off diag entries
  
  for i=0:39
    phi(i+(1),i+(1)) = 0.5*phi(i+(1),i+(1));
  end
     
  % pre-compute threshold for entering last loop
  % Equation (60)
  C3_sum = 0;
  C3_cnt = 0;
  max3   = 0;
  avg3   = 0;
  K3     = 0.4;  
  
  for i0=0:5:35
    for i1=1:5:36
      for i2=2:7:37
      
        C3_cnt += 1;
        C3      = dMag(i0+(1)) + dMag(i1+(1)) + dMag(i2+(1));
        C3_sum += C3;
        
        max3 = max(max3,C3);
      end
    end
  end  
  
  avg3 = C3_sum / C3_cnt;
  thr3 = avg3 + K3*(max3-avg3);
  
  % find best pulse positions, based on
  % maximizing Equation (53)

  maxRatio = 0;
  iPulse = [0,1,2,3];  % setup valid pulse positions

  % control number of times final loop entered
  L4_cnt = 0;  
  
  for i0=0:5:35        % pulse positions m0
    for i1=1:5:36      % pulse positions m1
      for i2=2:5:37    % pulse positions m2
        
        % correlation due to 3 pulses
        C = dMag(i0+(1)) + dMag(i1+(1)) + dMag(i2+(1));
     
        if( (C > thr3) && (L4_cnt < L4_max))
        
          L4_cnt += 1;
          for i3 = [3:5:38 4:5:39]  % pulse positions m3
            
            % compute value proportional to Ck^2 / E 
            E     = compute_E(phi,i0,i1,i2,i3);
            ratio = (C + dMag(i3+(1)))^2 / E;
            
            if( ratio > maxRatio)
              maxRatio = ratio;
              iPulse   = [i0 i1 i2 i3];
            end
            
          end % i3 loop      
        end % endif C > ...
        
      end  %i2 loop
    end % i1 loop
  end % i0 loop
  
  % based on best pulse positions, get the amplitude info
  sPulse = dSign(iPulse+(1));
  
  % output codebook vector
  c = zeros(1,40);
  c(iPulse+(1)) = sPulse;
 
  % account for pitch prefilter if delay < 40
  if( T < 40)
  
    N = 40 - T;
    
    % Equation (49)
    % add delayed version of h to last N samples
    c(end-N+1:end) = c(end-N+1:end) + B*c(0+(1):N-1+(1));
  end 
 
end

% represents Ek = ck'*PHI*ck, Energy using k'th possible codebook vector given
% by pulse positions (i0,i1,i2,i3)
function E = compute_E(phi,i0,i1,i2,i3)

  % Equation (59)
  E = phi(i0+(1),i0+(1)) +...
      phi(i1+(1),i1+(1)) + phi(i0+(1),i1+(1)) +... 
      phi(i2+(1),i2+(1)) + phi(i0+(1),i2+(1)) + phi(i1+(1),i2+(1)) +...
      phi(i3+(1),i3+(1)) + phi(i0+(1),i3+(1)) + phi(i1+(1),i3+(1)) + phi(i2+(1),i3+(1));
      
end

% correlation between x(i) and y(i-n);
function z = compute_correlation(x,y)

  assert(length(x) == length(y));
  L = length(x);
  z = zeros(1,L);
  
  for n=0:(L-1)
    sLen = L - n;
    s1 = x(end-sLen+1:end);     % last sLen samples
    s2 = y(0+(1):sLen-1+(1));   % first sLen samples
    z(n+(1)) = s1*s2'; 
  end
    
end