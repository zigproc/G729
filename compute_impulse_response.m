function h = compute_impulse_response(a,gamma1,gamma2,a_hat)
  % Compute imputse response
  
  assert(length(a) == 10);
  assert(length(a_hat) == 10);  % full impulse reponse is [1 a], and [1 a_hat]
  
  h = zeros(1,40);
  
  a1 = a.*(gamma1.^(1:10));
  a2 = fliplr(a.*(gamma2.^(1:10)));
  a_hat = fliplr(a_hat);
  
  % filter signal consisting of coefficients of A(z/gamma1) 
  % extended by zeros thru filters 1/A_hat(z) and 1/A(z/gamma2)
  
  x = zeros(1,40);
  x(1:11) = [1 a1];
  
  % [out(n-10) out(n-9) .. out(-1)]
  % TBD: is zero-state proper?
  state_fb = zeros(1,10);
  
  % filter thru 1/A_hat(z)
  for i = 1:length(x)
    sum = x(i);
    sum -= a_hat* state_fb';
    
    h(i) = sum;

    state_fb = [state_fb(2:10) h(i)];
  end
  
  % filter in-place thru 1/A(z/gamma2)
  state_fb = zeros(1,10);
  
  for i = 1:length(h)
    sum = h(i);
    sum -= a2 * state_fb';
    
    h(i) = sum;

    state_fb = [state_fb(2:10) h(i)];
  end  
  

end