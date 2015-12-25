function c = compute_fixed_codebook_vector(iPulse,sPulse,T,gp_prev)

   c = zeros(1,40);
   c(iPulse+(1)) = sPulse;
    
   if( T < 40)
      B = max(gp_prev,0.2);
      B = min(B,0.8);
    
      N = 40 - T;
    
      % Equation (48)
      % add delayed version of c to last N samples
      c(end-N+1:end) = c(end-N+1:end) + B*c(0+(1):N-1+(1));
   end   
end