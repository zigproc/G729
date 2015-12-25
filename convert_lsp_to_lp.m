function a = convert_lsp_to_lp(q)

  assert(length(q) == 10);
  
  a  = zeros(1,11);
  q1 = q(1:2:end);
  q2 = q(2:2:end);
  
  % matlab indices 1:7
  % spec indices -1,0,1,..5
  f1 = zeros(1,7);
  f2 = zeros(1,7);
  
  % array index of the 0th term
  ZI = 2;

  f1(ZI) = 1.0;
  f2(ZI) = 1.0;
  
  for i=1:5  
    Q1_i = q1(i);
    Q2_i = q2(i);
    
    f1(ZI+i) = -2*Q1_i*f1(ZI+i-1) + 2*f1(ZI+i-2);
    f2(ZI+i) = -2*Q2_i*f2(ZI+i-1) + 2*f2(ZI+i-2);
    
    for j=(i-1):-1:1
      f1(ZI+j) = f1(ZI+j) - 2*Q1_i*f1(ZI+j-1) + f1(ZI+j-2);
      f2(ZI+j) = f2(ZI+j) - 2*Q2_i*f2(ZI+j-1) + f2(ZI+j-2);
    end
    
  end
  
  f1Prime = f1;
  f2Prime = f2;
  
  % Equation (25)
  for i=1:5
    f1Prime(ZI+i) = f1(ZI+i) + f1(ZI+i-1);
    f2Prime(ZI+i) = f2(ZI+i) - f2(ZI+i-1);
  end
  
  for i=1:5
    a(i+1) = 0.5*(f1Prime(ZI+i) + f2Prime(ZI+i));
  end
 
  for i=6:10
    a(i+1) = 0.5*(f1Prime(ZI+11-i) - f2Prime(ZI+11-i));
  end
  
  a(1) = 1.0;
end