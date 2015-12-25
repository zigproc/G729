function [w] = convert_lp_to_lsp(a,x_os)

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
  
  r1 = find_roots(f1,x_os);
  r2 = find_roots(f2,x_os);

  q = zeros(1,10);
  q(1:2:end) = r1;
  q(2:2:end) = r2;
  w = acos(q);
  % check that roots interleave
  assert(length(find(diff(w) < 0)) == 0);
  
end

% find roots using either f1 or f2, and oversampled grid
function r = find_roots(f,x_os)

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
  
    start_idx = 4*(idx(i)-1) + (1);
    stop_idx  = 4*(idx(i))   + (1);
    % finer grid for tracking the root

    for j=start_idx:stop_idx
      x = x_os(j);
      q_fine(j - start_idx + (1)) = C(x,f);
    end
    
    root_idx = find(diff(sign(q_fine)) ~= 0);
    
    % in sub interval, should have one root
    assert(length(root_idx) == 1);
        
    % on 5 point sub interval grid, find indices between zero crossing
    start_idx = start_idx + root_idx - (1);
    
    % easy approx, split interval.
    r(i) = 0.5*(x_os(start_idx) + x_os(start_idx+1));
    
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
    
  