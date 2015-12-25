function [gainIdx,gp_hat_opt,gc_hat_opt,U_hat_opt] = quantize_gains(c,...
                                                                    x,...
                                                                    y,...
                                                                    z,...
                                                                    gc_opt,...
                                                                    gp_opt,...
                                                                    U_hat_prev,...
                                                                    gainCodebook)
   % shortcut
   cb = gainCodebook;
   
   % need to minimize Equation (63)   
   Ey  = y*y';
   Ez  = z*z';
   Exy = x*y';
   Exz = x*z';
   Eyz = y*z';
   
   % mean energy of fixed codebook contribution
   % Equation (66)
   E     = 10*log10((1/40)*(c*c'));
   Epred = cb.pred * U_hat_prev';
   Ebar  = 30;  %TBD seems to assume a particular range of optimal gc values
   
   % initialize for optimization loop
   Eopt = Inf;
   gainIdx.iA = 0;
   gainIdx.iB = 0;
   U_hat_opt  = 0;
   
   % full codebook search (TBD reduced complexity)
   for iA=1:length(cb.gA)
     for iB=1:length(cb.gB)
     
        % build the candidate quantized ACB gain from the CB
        gp_hat_cand      = cb.gA(iA,1) + cb.gB(iB,1);
        % build the candidate quantized correction factor from the CB
        gamma_cand       = cb.gA(iA,2) + cb.gB(iB,2);
        U_hat_cand       = 20*log10(gamma_cand);
        
        % gc in log domain after accouting for prediction, mean subtract,
        % and correction factor
        gc_hat_cand_log  = Epred + Ebar - E + U_hat_cand;
        % convert gc to linear domain
        gc_hat_cand      = 10^(gc_hat_cand_log/20);
        
        % Equation (63)
        % don't need to include Ex = x*x'
        Ecand = (gp_hat_cand^2)*Ey  +...
                (gc_hat_cand^2)*Ez  -...
                (2*gp_hat_cand)*Exy -...
                (2*gc_hat_cand)*Exz +...
                (2*gp_hat_cand*gc_hat_cand)*Eyz;
                 
        % minimizing error
        if(Ecand < Eopt)
          
          gainIdx.iA = iA-(1);  % converting to indices for transmission
          gainIdx.iB = iB-(1);
          U_hat_opt  = U_hat_cand;
          Eopt       = Ecand;
          gp_hat_opt = gp_hat_cand;
          gc_hat_opt = gc_hat_cand;
          
        end       
     end   
   end
                                                  
end