function y = rearrange(x,J)
% rearrange candidate vectors.
% x is DxN, where D is dimension, and N number of vectors
% J is min distance

% section 3.2.4
  [D,N] = size(x);
  
  y = zeros(D,N);
  
  for i=1:N
    v = x(:,i);
    for k=2:D
      if(v(k-1) > v(k) - J)
        v(k-1) = (v(k) + v(k-1) - J)/2;
        v(k)   = (v(k) + v(k-1) + J)/2;
      end
    end
    
    y(:,i) = v;
end
