function [s_dec,decState,signals] = g729_decode( decState, pkt, consts, rom, cb)
  
    Nsf = consts.subframeSize;
    Nf  = 2*Nsf;                      % framesize
    M   = consts.lpcOrder;  
  
    signals = [];
    
    % if packet contains bitsream, decode, and replace pkt with params
    % extracted from bitstream.  Should match!!
    if( isfield(pkt,"bitstream"))
      pkt = g729_decode_bitstream(pkt.bitstream);
    end
    
    % interpolate LSP and convert to synthesis filter
    if( isfield(pkt,"lsfIdx"))  % check if LSF sent as quantized codebook entries
      [l_hat,wi_hat] = decode_lsf(pkt.lsfIdx, cb.lsf, decState.li_hat_prev);
      decState.li_hat_prev = [decState.li_hat_prev(:,2:end) l_hat]; 
      wi_hat = wi_hat';
    else
      wi_hat = pkt.wi_hat;
    end  
      
    qi_hat = cos(wi_hat);
    
    lsp_hat_subframe{1} = 0.5*(qi_hat + decState.qi_hat_prev);
    lsp_hat_subframe{2} = qi_hat;
    decState.qi_hat_prev = qi_hat;
    
    a_hat{1} = convert_lsp_to_lp(lsp_hat_subframe{1});
    a_hat{2} = convert_lsp_to_lp(lsp_hat_subframe{2});

    for subframeIdx=1:2

      % Decode adaptive codebook vector
      T_int    = pkt.sf{subframeIdx}.T_int;
      T_frac   = pkt.sf{subframeIdx}.T_frac;
      ZI_U     = decState.ZI_U;
      ZI_S_HAT = decState.ZI_S_HAT;
    
      % Decode ACB
      v = compute_adaptive_codebook_vector(decState.u,...
                                           ZI_U,...
                                           T_int,...
                                           T_frac,...
                                           rom.interp_u);
                                       
      % Decode fixed codebook vector
      K = 2^13;  % tmp
      c = compute_fixed_codebook_vector(pkt.sf{subframeIdx}.iPulse,...
                                        pkt.sf{subframeIdx}.sPulse,...
                                        T_int,...
                                        decState.gp_hat_prev);
      
      % check if quantized gain indices available.  decode if so.
      if( isfield(pkt.sf{subframeIdx},"gainIdx"))
        [gp_hat,gc_hat,U_hat] = decode_gains(pkt.sf{subframeIdx}.gainIdx.iA + (1), % add 1 since iA is tx'd param, 0 indexed
                                             pkt.sf{subframeIdx}.gainIdx.iB + (1),
                                             c,
                                             decState.U_hat_prev,
                                             cb.gain);
                                             
         decState.U_hat_prev = [decState.U_hat_prev(2:end) U_hat];
      else
         gp_hat = pkt.sf{subframeIdx}.gp_hat;
         gc_hat = pkt.sf{subframeIdx}.gc_hat;
      end 

 
      % Compute excitation for this subframe
      u = (gp_hat)*v  + (gc_hat)*c/K;
      
      decState.u(ZI_U:ZI_U+Nsf-1) = u;
    
      % compute reconstructed speech
      sf_state.ahat.fb = decState.s_hat(decState.ZI_S_HAT-10:ZI_S_HAT-1);
      [s_hat,~] = compute_synthetic_speech(u,...
                                           sf_state,...
                                           a_hat{subframeIdx}(2:11));
                                           
      decState.s_hat(ZI_S_HAT:ZI_S_HAT+Nsf-1) = s_hat;
      
      % postfilter the speech                                     
      [s_post,res,res_pf,decState] = dec_postfilter( decState,...
                                                     a_hat{subframeIdx}(2:11),...
                                                     Nsf,
                                                     T_int,...
                                                     rom.gammaN,...
                                                     rom.gammaD,...
                                                     rom.gammaP,...
                                                     rom.interp_rhat); 
      
      % highpass filter the speech 
      [s_dec, decState.postproc.h2] = preproc(s_post,...
                                              decState.postproc.h2,...
                                              rom.h2.b,...
                                              rom.h2.a); 
      % restore signal level
      s_dec = 2*s_dec;     
            
      % update gp_hat_prev
      decState.gp_hat_prev = gp_hat;
      
      % shift down states for next subframe
      Nkeep = length(decState.rhat) - Nsf;
      decState.rhat = [decState.rhat(end-Nkeep+1:end) zeros(1,Nsf)];
  
      Nkeep = length(decState.s_hat) - Nsf;
      decState.s_hat = [decState.s_hat(end-Nkeep+1:end) zeros(1,Nsf)];
      
      Nkeep = length(decState.u) - Nsf;
      decState.u = [decState.u(end-Nkeep+1:end) zeros(1,Nsf)];
           
      % save interesting signals
      
      r1 = (subframeIdx-1)*Nsf + (1);
      r2 = r1 + Nsf - 1;
      
      signals.u(r1:r2)      = u;
      signals.s_hat(r1:r2)  = s_hat;
      signals.s_post(r1:r2) = s_post;
      signals.s_dec(r1:r2)  = s_dec;
      signals.residual_dec(r1:r2) = res;
      signals.residual_pf_dec(r1:r2) = res_pf; 
    
    end
end