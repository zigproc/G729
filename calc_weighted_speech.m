function [sw,ws_state] = calc_weighted_speech(s,ZI_S,N,ws_state,a,gamma1,gamma2)

  % ws_state.ff = [s(n-10) s(n-9) .. s(n-1)]
  % ws_state.fb = [sw(n-10) sw(n-9) .. sw(n-1)]
  
  % Perceptual weighting filter W(z) = A(z/gamma1) / A(z/gamma2)
  a_ff = fliplr(a .* (gamma1 .^ (1:10)));
  a_fb = fliplr(a .* (gamma2 .^ (1:10)));
    
  assert(ZI_S+N-1 <= length(s));
  
  sw = zeros(1,N);
  
  for i=1:N
  
     sample = s(ZI_S+i-1);
     % Equation (33)
     sum  = sample;
     sum += a_ff * ws_state.ff';
     sum -= a_fb * ws_state.fb';
     
     sw(i) = sum;
     
     % update history
     ws_state.ff = [ws_state.ff(2:10) sample];
     ws_state.fb = [ws_state.fb(2:10) sw(i)];
  
  end
  

end