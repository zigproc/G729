function [x,residual] = compute_target_signal(curr_subframe,r_state,a_hat,gamma1,gamma2,x_state)

  % Section 3.6
  % target signal computed by filtering the LP residual signal  
  % r(n) through the combination of the synthesis filter (1/A_hat(z)) and
  % the weighting filter A(z/gamma1) / A(z/gamma2)
  
  % r_state = [s(n-10) s(n-9) ... s(n-1)]
  
  assert(length(curr_subframe) == 40);
  assert(length(a_hat) == 10);
  assert(length(r_state.ahat.ff) == 10);
  assert(length(x_state.w.ff) == 10);
  assert(length(x_state.w.fb) == 10);
  assert(length(x_state.ahat.fb) == 10);
  
  residual = zeros(1,40);
  x = zeros(1,40);  % target signal
  
  a1 = fliplr(a_hat.*(gamma1.^(1:10)));  
  a2 = fliplr(a_hat.*(gamma2.^(1:10)));
  a_hat = fliplr(a_hat);
  
  % compute residual signal
  % equation (36)
  for i=1:length(curr_subframe)
  
    sample = curr_subframe(i);
    sum  = sample;
    sum += a_hat * r_state.ahat.ff';  
    
    % output
    residual(i) = sum;
    
    % update FF state (past input values s(n))
    r_state.ahat.ff = [r_state.ahat.ff(2:10) sample];
    
  end
  
  % filter residual signal through synthesis filter
  
  for i=1:length(residual)
    sum  = residual(i);
    sum -= a_hat * x_state.ahat.fb';
    
    x(i) = sum;
    % update FB state (past output values x
    x_state.ahat.fb = [x_state.ahat.fb(2:10) sum];
  end
  
  
  % filter x through W(z)
  for i=1:length(x)

    sample = x(i);
    sum = sample;
    sum -= a2 * x_state.w.fb';
    sum += a1 * x_state.w.ff';
    x(i) = sum;
    
    x_state.w.fb = [x_state.w.fb(2:10) sum];
    x_state.w.ff = [x_state.w.ff(2:10) sample];
  end
  
end