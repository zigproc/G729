function [q] = convert_lp_to_lsp_bisect(a,x_os,num_bisections)

  % +(1) generally indicates when a given formula in the spec is 0-indexed,
  % so this is used to convert to matlab indexing.  Other formulas in the spec
  % that are 1-indexed are used directly in Matlab.  Need to be mindful when 
  % converting equations for C 

  assert(length(a) == 11);

  f1 = zeros(1,6);
  f2 = zeros(1,6);
  f1(1) = 1.0;
  f2(1) = 1.0;

  % Equation (15)
  % f1(i+1) = a(i+1) + a(10-i) - f1(i)  
  % f2(i+1) = a(i+1) - a(10-i) - f2(i)  

  % since we grab a1--a5, and a10--a6 based on spec indexing, this
  % alters accesses thru af/ab. 
  af = a(2:6);
  ab = fliplr(a(7:11));
  
  for i = 0:4
    f1(i+1+(1)) = af(i+(1)) + ab(i+(1)) - f1(i+(1));
    f2(i+1+(1)) = af(i+(1)) - ab(i+(1)) + f2(i+(1));
  end
  
  r1 = find_roots(f1,x_os,num_bisections);
  r2 = find_roots(f2,x_os,num_bisections);

  q = zeros(1,10);
  q(1:2:end) = r1;
  q(2:2:end) = r2;
  
  w = acos(q);
  % check that roots interleave
  assert(length(find(diff(w) < 0)) == 0);
  
end

% find roots using either f1 or f2, and oversampled grid
function r = find_roots(f,x_os,num_bisections)

 
  % currently evaluate the function at 60 points 
  q = zeros(1,60);
  
  % main points from oversampled grid
  x_coarse = x_os(1:4:end);
  assert(length(x_coarse) == 60); 
  
  % evaluate at 60 points
  % use either f1 or f2
  for i = 1:length(x_coarse)
    x = x_coarse(i);
    q(i) = C(x,f);
  end
  
  % find zero crossings
  idx = find(diff(sign(q)) ~= 0);
  
  % expect 5 roots
  assert(length(idx) == 5);
  
  r = zeros(1,5);
  
  q_fine = zeros(1,5);
  
  % TBD: need to use bisection instead
  for i=1:length(idx)
  
    xLo = x_coarse(idx(i));    
    xHi = x_coarse(idx(i)+1);
    
    yLo = C(xLo,f);
    yHi = C(xHi,f);
    
    % keep subdividing in half where root is located
    for b=1:num_bisections
            
      xMid = 0.5*(xLo+xHi);
      yMid = C(xMid,f);
      
      % midpoint is the root, so exit bisection loop
      if( yMid == 0)
        break;
      end
      
      % check which half the root is in, and
      % update boundaries accordingly
      if( sign(yLo) ~= sign(yMid))
        xHi = xMid;
        yHi = yMid;
      elseif( sign(yMid) ~= sign(yHi))
        xLo = xMid;
        yLo = yMid;
      else
        assert(false,'Error in bisection');
      end
      
    end
    
    if( yMid == 0)
      r(i) = xMid;
    else
      % linearly interpolate intercept based on boundary points
      invm = (xHi - xLo) / (yHi - yLo);
      r(i) = xLo - yLo*invm;
    end  
    
  end
  
end  
% Equation (17)
% Evaluate C(x) = T5(x) + f(1)T4(x) + f(2)T3(x) + f(3)T2(x) + f(4)T1(x) + f(5)/2
% evaluated using recursive relationship given in spec
function y = C(x,f)
  
  % Assume f length 6.  Spec uses 0-indexing for f, but 1-indexing for b

  b = zeros(1,6);
  b(6) = 0;
  b(5) = 1.0;
  
  for k = 4:-1:1
    b(k) = 2*x*b(k+1) - b(k+2) + f(5-k+(1));
  end
  
  y = x*b(1) - b(2) + 0.5*f(5+(1));
  
end
    
  