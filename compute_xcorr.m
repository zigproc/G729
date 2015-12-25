function r = compute_xcorr(s1,s2,nLags)
  
  % cross correlation, <s(i),s(i-k)>

  assert(length(s1) == length(s2));
  assert(length(s1) >= nLags);
  
  N = length(s1);

  r = zeros(1,nLags);
  
  for k = 0:(nLags-1)
    
    % Equation (5)
    sLen  = N - k;
    s1_x  = s1(k+(1):(N-1+(1)));        % last sLen samples
    s2_x  = s2(0+(1):sLen-1+(1));       % first sLen samples
    
    r(k+1) = s1_x*s2_x';
   end

end
