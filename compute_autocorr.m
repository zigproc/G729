function r = compute_autocorr( s, w_lp, w_lag, lag0_Min, noiseFloor)
% compute autocorrelation of speech segment "s" using 
% window function w_lp, modified by bandwidth expansion

Nsamples = length(s);

% window and input signal need to be same size
assert(length(s) == length(w_lp))

% w_lag modifies all but first autocorrelation lag
% total length of r is therefore length(w_lag)+1
Nlags = length(w_lag) + 1;
r = zeros(1,Nlags);

% G729 Equation (4)
s_prime = w_lp .* s;

 
r = compute_xcorr(s_prime,s_prime,Nlags);

# lower boundary for lag 0
r(1) = max( r(1), lag0_Min);

# Noise floor and bandwidth expansion
# Equation (7)
r(1) = noiseFloor*r(1);
r(2:end) = r(2:end) .* w_lag;

 
  

