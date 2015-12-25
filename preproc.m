function [s,state] = preproc(input,state,b,a)

  assert(length(input) == 80 || length(input) == 40);
  assert(length(state.fb) == 2);
  assert(length(state.ff) == 2);
  assert(length(b) == 3);
  assert(length(a) == 2);
  
  a = fliplr(a); 
  b = fliplr(b);
  
  s = zeros(1,length(input));
  
  for n=0:(length(input)-1)
    sample = input(n+(1));
    
    % augmented state with current sample
    ff = [state.ff sample];
    % 2nd order IIR
    sum  = ff * b';
    sum -= state.fb * a';
    
    % output, and update state
    s(n+(1)) = sum;
    state.ff = [state.ff(2:end) sample];
    state.fb = [state.fb(2:end) sum]; 
  end

end