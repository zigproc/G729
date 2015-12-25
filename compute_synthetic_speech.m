function [s_hat,syn_state] = compute_synthetic_speech(u,syn_state,a_hat)

  assert(length(u) == 40);
  assert(length(a_hat) == 10);
  assert(length(syn_state.ahat.fb) == 10);
  
  a_hat = fliplr(a_hat);
  
  for n=0:39
    sample = u(n+(1));
    sum    = sample;
    sum   -= a_hat * syn_state.ahat.fb';
    
    s_hat(n+(1)) = sum;
    
    syn_state.ahat.fb = [syn_state.ahat.fb(2:end) sum];
    
end
