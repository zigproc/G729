function [gp_hat,gc_hat,U_hat] = decode_gains(iA,iB,c,U_hat_prev,gainCodebook)
   % shortcut
   cb = gainCodebook;
   
   % mean energy of fixed codebook contribution
   % Equation (66)
   E     = 10*log10((1/40)*(c*c'));
   Epred = cb.pred * U_hat_prev';
   Ebar  = 30;  %TBD seems to assume a particular range of optimal gc values
   
   % build the candidate quantized ACB gain from the CB
   gp_hat = cb.gA(iA,1) + cb.gB(iB,1);
   % build the candidate quantized correction factor from the CB
   gamma  = cb.gA(iA,2) + cb.gB(iB,2);
   U_hat  = 20*log10(gamma);
  
   % gc in log domain after accouting for prediction, mean subtract,
   % and correction factor
   gc_hat_log  = Epred + Ebar - E + U_hat;
   % convert gc to linear domain
   gc_hat      = 10^(gc_hat_log/20);
                                                  
end