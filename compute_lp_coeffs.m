function [a,k, E] = compute_lp_coeffs( r)

# Solve Equation (8)
N = length(r);
a = zeros(1,N);
a_tmp = zeros(1,N);
k = zeros(1,N-1);

E = r(1);
a(1) = 1;


for i = 1:(N-1)

  a_tmp = a;
  
  s1 = a(1:i);     % a coeffs 1 to i in forward order
  s2 = r(2:(i+1)); 
  s2 = fliplr(s2); % lag 1 to i in reverse order
  
  % reflection coeff
  k(i) = (-s1*s2') / E;
  
  % newest coeff
  a(i+1) = k(i);
  
  % modify inner coeffs
  for j = 1:(i-1)
    a(j+1) = a(j+1) + k(i)*a_tmp(i-j+1);
  end
  

  E = (1- k(i)^2)*E;
  
end



