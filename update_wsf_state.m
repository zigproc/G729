function wsf_state = update_wsf_state(wsf_state,input,ew,a_hat)
 
  % filter input thru synthesis filter
  assert(length(wsf_state.ahat.fb) == length(a_hat));
  assert(length(wsf_state.w.ff) == length(a_hat));
  assert(length(wsf_state.w.fb) == length(a_hat));
  assert(length(input) == 40);
  assert(length(ew) == 40);
  
  a_hat = fliplr(a_hat);
  e     = zeros(1,length(input));
  
  for i=1:length(input)
    sum  = input(i);
    sum -= a_hat * wsf_state.ahat.fb';
    
    e(i) = sum;
    % update FB state (past output values x
    wsf_state.ahat.fb = [wsf_state.ahat.fb(2:10) sum];
  end
  
  % Use e(n) and ew(n) as if e(n) went thru weighted synthesis filter
  Nkeep = length(a_hat);
  wsf_state.w.ff = e(end-Nkeep+1:end);
  wsf_state.w.fb = ew(end-Nkeep+1:end);
    
end