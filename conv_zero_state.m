function [y] = conv_zero_state(h,x)

  assert(length(h) == length(x));
  % convolve signal x with filter h
  % zero-state convolution  y(n) = sum k=0 to n { h(k)*x(n-k) }
  
  y = conv(h,x);
  y = y(1:length(x));
  
end