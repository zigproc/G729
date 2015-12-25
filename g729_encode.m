% clear 
clear all;
close all;

DECODE_ENABLE = 1;
LSF_QUANTIZATION_ENABLE = 1;
GAIN_QUANTIZATION_ENABLE = 1;
DEBUG_ENABLE = 1;

% signal to test
filename = 'OSR_us_000_0030_8k.wav';

% init the G729 ROM tables
init_rom;

% init the G729 encoder state
init_ram_enc;

if( DECODE_ENABLE)
  init_ram_dec;
end

if(LSF_QUANTIZATION_ENABLE || GAIN_QUANTIZATION_ENABLE || DECODE_ENABLE)
  init_codebook_rom;
end

% signal to use
[yFile,fs,bps]=wavread(filename);
yFile  = yFile';
yi = 1;

% grab some interesting segment
yStart = 46880;
yLen   = 16000;
yFile = yFile(yStart:yStart+yLen-1);

nFrames = floor(length(yFile) / g729.consts.frameSize);
frameIdx = 0;

% save some frame by frame stats
frameStats = [];
frameStats.energy = zeros(1,nFrames);
frameStats.k      = zeros(10,nFrames);
frameStats.T_op   = zeros(1,nFrames);
frameStats.wi     = zeros(10,nFrames);
frameStats.r_corr = zeros(11,nFrames);

% shortcuts
Nsf = g729.consts.subframeSize;
Nf  = 2*Nsf;                      % framesize
M   = g729.consts.lpcOrder;

if DEBUG_ENABLE
  % signals to save
  signals.enc.residual = zeros(1,Nf*nFrames); % residual
  signals.enc.x        = zeros(1,Nf*nFrames); % target vector prior to ACB search
  signals.enc.u        = zeros(1,Nf*nFrames); % excitation
  signals.enc.c        = zeros(1,Nf*nFrames); % FCB vector (unscaled)
  signals.enc.v        = zeros(1,Nf*nFrames); % ACB vector (unscaled)
  signals.enc.gp_hat   = zeros(1,Nf*nFrames); % ACB gain
  signals.enc.gc_hat   = zeros(1,Nf*nFrames); % FCB gain
  signals.enc.gp_opt   = zeros(1,Nf*nFrames); % ACB gain
  signals.enc.gc_opt   = zeros(1,Nf*nFrames); % FCB gain
  signals.enc.s        = zeros(1,Nf*nFrames); % preprocessed speech
  signals.enc.sw       = zeros(1,Nf*nFrames); % weighted speech
  signals.enc.s_hat    = zeros(1,Nf*nFrames); % synthesized speech
  signals.enc.ew       = zeros(1,Nf*nFrames); % weighted error signal
  signals.enc.h        = zeros(1,Nf*nFrames); % impulse response
  signals.enc.mse_opt  = zeros(1,Nf*nFrames); % MSE when using optimal gains
  signals.enc.mse_quant= zeros(1,Nf*nFrames); % MSE when using quantized gains
  signals.enc.U_hat_opt= zeros(1,Nf*nFrames); % quantized FCB gain adj gamma (log domain)

  signals.dec.u      = zeros(1,Nf*nFrames);
  signals.dec.s_hat  = zeros(1,Nf*nFrames);
  signals.dec.s_post = zeros(1,Nf*nFrames);
  signals.dec.out    = zeros(1,Nf*nFrames);
  signals.dec.residual    = zeros(1,Nf*nFrames);
  signals.dec.residual_pf = zeros(1,Nf*nFrames);
end


% packets (TBD real quantization of LSF and gains)
pkt = cell(1,nFrames);


while true
    
  % get frame of data
  % break if can't get a full frame
  try
    input = yFile(yi:yi+g729.consts.frameSize-1);
    frameIdx += 1;
    yi += g729.consts.frameSize;
  catch
    break;
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %                       %
  % FRAME PROCESSING      %
  %                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Preprocess the input signal
  [s_curr, g729.enc.state.preproc.h1] = preproc(input,...
                                                g729.enc.state.preproc.h1,...
                                                g729.enc.h1.b,...
                                                g729.enc.h1.a);
  
  % shift in new frame of data to window
  Nkeep = length(g729.enc.state.speechBuf) - length(s_curr);
  g729.enc.state.speechBuf = [g729.enc.state.speechBuf(end-Nkeep+1:end) s_curr];
  
  % shift weighted speech down frame
  Nkeep = length(g729.enc.state.weightedSpeechBuf) - g729.consts.frameSize;
  g729.enc.state.weightedSpeechBuf = [g729.enc.state.weightedSpeechBuf(end-Nkeep+1:end) zeros(1,g729.consts.frameSize)];
  
  
  % compute autocorrelation with -40 dB noise floor
  % Section 3.2.1
  r = compute_autocorr(g729.enc.state.speechBuf,g729.enc.w_lp,g729.enc.w_lag,0.0,1.0001);

  % run levinson durbin algorithm
  % Section 3.2.2
  [a_lp,k,E] = compute_lp_coeffs( r);
  
  % convert LP to LSF
  % Section 3.2.3
  qi = convert_lp_to_lsp_bisect(a_lp,g729.enc.lsp_qi,4);   % unquantized LSP
  wi = acos(qi);      % unquantized LSF 
  
  % fake quantization of LSF
  % section 3.2.4
if LSF_QUANTIZATION_ENABLE
  [l_hat,wi_hat,idx] = quantize_lsf(wi,...
                                    g729.enc.state.li_hat_prev,...
                                    g729.rom.cb.lsf);
                                   
  % shift in latest lhat
  g729.enc.state.li_hat_prev = [g729.enc.state.li_hat_prev(:,2:end) l_hat]; 
  wi_hat = wi_hat';  
else 
  wi_hat = wi;
end

  qi_hat = cos(wi_hat);
  
  % interpolation.  Section 3.2.5
  % interpolate the unquantized LSP
  % Equation (24)
  lsp_subframe{1} = 0.5*(qi + g729.enc.state.qi_prev);
  lsp_subframe{2} = qi;
  g729.enc.state.qi_prev = qi;
 
  % interpolate the quantized LSP
  lsp_hat_subframe{1} = 0.5*(qi_hat + g729.enc.state.qi_hat_prev);
  lsp_hat_subframe{2} = qi_hat;
  g729.enc.state.qi_hat_prev = qi_hat; 
  
  % unquantized, interpolated LP coeffs
  a{1} = convert_lsp_to_lp(lsp_subframe{1});
  a{2} = a_lp;
  
  % LSP to LP conversion
  % using interpolated values of quantized LSP
  % Section 3.2.6
  a_hat{1} = convert_lsp_to_lp(lsp_hat_subframe{1});
  a_hat{2} = convert_lsp_to_lp(lsp_hat_subframe{2});
  
  % Equation (28)
  LAR = (1 + k(1:2)) ./ (1 - k(1:2));
  LAR = log10(LAR);  
  
  % interpolate the LAR
  LAR_subframe{1} = 0.5*(LAR + g729.enc.state.LAR_prev);
  LAR_subframe{2} = LAR;
  g729.enc.state.LAR_prev = LAR;
  
  % save pkt data (frame based)
  
if LSF_QUANTIZATION_ENABLE
  pkt{frameIdx}.lsfIdx = idx;
else
  pkt{frameIdx}.wi_hat = wi_hat; 
end 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %                       %
  % SUBFRAME PROCESSING   %
  % (Part 1)              %
  %                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%  

  % shortcuts
  Nsf = g729.consts.subframeSize;
  s = g729.enc.state.speechBuf;

  gamma = {};
  
  for subframeIdx=1:2

    ws = g729.enc.state.weightedSpeechBuf;   % sw(n) in spec!
  
    % indices into speech and weighted speech, current subframe
    ZI_S_SF  = g729.enc.state.ZI_S  + (subframeIdx-1)*Nsf;    % zero index of current subframe in speechBuf
    ZI_WS_SF = g729.enc.state.ZI_WS + (subframeIdx-1)*Nsf;    % zero index of current subframe in wgt speechBuf
    
    % Perceptual weighting
    % Section 3.3
    % update value of flat
    [gamma1,gamma2,g729.enc.state.flat_prev] = compute_gamma(LAR_subframe{subframeIdx},...
                                                             g729.enc.state.flat_prev,...
                                                             wi);
    
    % save gamma values, which define W(z) each subframe, along with subframe a_hat    
    gamma{subframeIdx}.gamma1 = gamma1;
    gamma{subframeIdx}.gamma2 = gamma2;    
  
    % history buffers for W(z) filter
    ws_state.ff = s(ZI_S_SF-10:ZI_S_SF-1);        % 10 samples, from prev subframes
    ws_state.fb = ws(ZI_WS_SF-10:ZI_WS_SF-1);     % 10 samples from prev weighted speech   
  
    % calculate weighted speech for current subframe
    [ws_curr,ws_state] = calc_weighted_speech(g729.enc.state.speechBuf,...
                                              ZI_S_SF,...
                                              Nsf,...
                                              ws_state,...
                                              a{subframeIdx}(2:11),...
                                              gamma1,gamma2);
    
    % weighted speech n=0...39 in current subframe
    g729.enc.state.weightedSpeechBuf(ZI_WS_SF:ZI_WS_SF+Nsf-1) = ws_curr;


if DEBUG_ENABLE      
    % DONE!    
    % save interesting signals                                     
    r1 = (frameIdx-1)*Nf + (subframeIdx-1)*Nsf + 1;
    r2 = r1 + Nsf - 1;
  
    signals.s(r1:r2)  =  g729.enc.state.speechBuf(ZI_S_SF:ZI_S_SF+Nsf-1);
    signals.sw(r1:r2) =  ws_curr;  
end
    
  end
  
  % Open Loop pitch analysis
  % section 3.4
  T_op = open_loop_pitch_analysis(g729.enc.state.weightedSpeechBuf,...
                                  g729.enc.state.ZI_WS);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %                       %
  % SUBFRAME PROCESSING   %
  % (Part 2)              %
  %                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%  

  s = g729.enc.state.speechBuf;
  
  % center point of ACB search
  % for first subframe it is the open loop pitch period
  T_acb_center = T_op;

  % max number of innermost loop of codebook search (pulse i3), across both subframes
  L4_max = 180;
  
  for subframeIdx=1:2

    % indices into speech and weighted speech, current subframe
    ZI_S_SF  = g729.enc.state.ZI_S  + (subframeIdx-1)*Nsf;    % zero index of current subframe in speechBuf
    ZI_U     = g729.enc.state.ZI_U;
    
    gamma1 = gamma{subframeIdx}.gamma1;
    gamma2 = gamma{subframeIdx}.gamma2;
    
    % compute impulse response
    % Section 3.5
    % use gamma1/gamma2 on subframe basis, as well as quantized/interpolated synthesis filter (LP)
    h = compute_impulse_response(a{subframeIdx}(2:11),...
                                 gamma1,gamma2,...
                                 a_hat{subframeIdx}(2:11)); 
   
    % Compute target signal
    % section 3.6
    curr_subframe = s(ZI_S_SF:ZI_S_SF+Nsf-1);   % 40 samples, current subframes
    r_state.ahat.ff = s(ZI_S_SF-10:ZI_S_SF-1);  % 10 samples, from prev subframes


    [x,residual] = compute_target_signal(s(ZI_S_SF:ZI_S_SF+Nsf-1),...
                                         r_state,...
                                         a_hat{subframeIdx}(2:11),...
                                         gamma1,gamma2,...
                                         g729.enc.state.wsf);
                                     
    % Section 3.7
    %copy LP residual into u for adaptive codebook search
    g729.enc.state.u(ZI_U:ZI_U+Nsf-1) = residual;
    
    [T_int,T_frac] = adaptive_codebook_search(T_acb_center,...
                                              subframeIdx,...  
                                              x,...
                                              g729.enc.state.u,...
                                              ZI_U,...
                                              h,...
                                              g729.enc.acb.interp_r);
    
    % starting point for ACB search in 2nd subframe is integer part for first subframe
    T_acb_center = T_int;                             
      
    % Section 3.7.1
    v = compute_adaptive_codebook_vector(g729.enc.state.u,...
                                         ZI_U,...
                                         T_int,...
                                         T_frac,...
                                         g729.enc.acb.interp_u);

    % Section 3.7.3 
    [gp,y] = compute_adaptive_codebook_gain(v,x,h);

    % Modify target signal
    % Equation (50)
    x_fcb = x - gp*y;

    % Fixed codebook search
    % Section 3.8
    [iPulse,sPulse,c,L4_cnt] = fixed_codebook_search(h,...
                                                     T_int,...
                                                     x_fcb,...
                                                     g729.enc.state.gp_hat_prev,...
                                                     L4_max); 
                                                   
    L4_max -= L4_cnt;
  
    assert(L4_max >= 0);  
                                             
    % TBD gain quantization  
    % Section 3.9
    K = 2^13;  % tmp to make sure gc is in proper range for gain quantization (hacky?!)
    z = conv_zero_state(c/K,h);
    
    [gp_opt,gc_opt] = compute_optimal_codebook_gains(x,x_fcb,y,z);

    % check    
    Eopt = gain_quantization_calc_MSE(x,y,z,gp_opt,gc_opt);
   
if GAIN_QUANTIZATION_ENABLE   
    [gainIdx,gp_hat,gc_hat,U_hat_opt] = quantize_gains(c,...
                                                       x,...
                                                       y,...
                                                       z,...
                                                       gc_opt,...
                                                       gp_opt,...
                                                       g729.enc.state.U_hat_prev,...
                                                       g729.rom.cb.gain);
                                                       
    % shift in newest quantized gc correction factor                                                   
    g729.enc.state.U_hat_prev = [g729.enc.state.U_hat_prev(2:end) U_hat_opt];
    
    Equant = gain_quantization_calc_MSE(x,y,z,gp_hat,gc_hat);
else    
    gp_hat = gp_opt;  % fake quantize
    gc_hat = gc_opt;
end    
    % update previous quantized ACB gain
    g729.enc.state.gp_hat_prev = gp_hat;

    % Memory Update, excitation buffer
    % Section 3.10
    utmp = gp_hat*v + gc_hat*(c/K);
    g729.enc.state.u(ZI_U:ZI_U+Nsf-1) = utmp;    

    % r(n) - u(n) in current subframe
    % this signal is filtered by synthesis filter update the memory 
    input = residual - utmp;
    % weighted error signal
    ew    = x - gp_hat*y - gc_hat*z;  
    
    % update the states of the weighted synthesis filter
    g729.enc.state.wsf = update_wsf_state(g729.enc.state.wsf,...
                                          input,...
                                          ew,...
                                          a_hat{subframeIdx}(2:11));

    % generate synthesized speech for testing
    [s_hat,g729.enc.state.sf] = compute_synthetic_speech(utmp,...
                                                         g729.enc.state.sf,...
                                                         a_hat{subframeIdx}(2:11));
    
    % shift excitation buffer down Nsf samples     
    Nkeep = length(g729.enc.state.u) - Nsf;
    g729.enc.state.u = [g729.enc.state.u(end-Nkeep+1:end) zeros(1,Nsf)];
    
    % save packet info for this subframe
    pkt{frameIdx}.sf{subframeIdx}.iPulse = iPulse;
    pkt{frameIdx}.sf{subframeIdx}.sPulse = sPulse;
    pkt{frameIdx}.sf{subframeIdx}.gp_hat = gp_hat;
    pkt{frameIdx}.sf{subframeIdx}.gc_hat = gc_hat;
    pkt{frameIdx}.sf{subframeIdx}.gp_opt = gp_opt;
    pkt{frameIdx}.sf{subframeIdx}.gc_opt = gc_opt;    
    pkt{frameIdx}.sf{subframeIdx}.T_int  = T_int;
    pkt{frameIdx}.sf{subframeIdx}.T_frac = T_frac;
    
if GAIN_QUANTIZATION_ENABLE   
    pkt{frameIdx}.sf{subframeIdx}.gainIdx = gainIdx;
end 
 
if DEBUG_ENABLE
    % DONE!    
    % save interesting signals                                     
    r1 = (frameIdx-1)*Nf + (subframeIdx-1)*Nsf + 1;
    r2 = r1 + Nsf - 1;    

    % save residual
    signals.enc.residual(r1:r2) = residual;
    % save target signal
    signals.enc.x(r1:r2) = x; 
    % save excitation signal
    signals.enc.u(r1:r2) = utmp;
    % save unscaled FCB signal
    signals.enc.c(r1:r2) = c; 
    % save FCB gain
    signals.enc.gc_hat(r1:r2) = gc_hat;
    % save FCB gain
    signals.enc.gc_opt(r1:r2) = gc_opt;
    % save unscaled ACB signal
    signals.enc.v(r1:r2) = v;
    % save ACB gain
    signals.enc.gp_hat(r1:r2) = gp_hat;
    % save ACB gain
    signals.enc.gp_opt(r1:r2) = gp_opt;
    % save synthetic speech
    signals.enc.s_hat(r1:r2) = s_hat;
    % save weighted error
    signals.enc.ew(r1:r2) = ew;
    % saved weighted impulse response
    signals.enc.h(r1:r2) = h;
    % save MSE
    signals.enc.mse_opt(r1:r2) = Eopt;
end    
    
    
if GAIN_QUANTIZATION_ENABLE 
    % save MSE
    signals.enc.mse_quant(r1:r2) = Equant;
    % save correction 
    signals.enc.U_hat_opt(r1:r2) = U_hat_opt;
end    

  
  end  % end subframe loop

  % format the bitstream  
if GAIN_QUANTIZATION_ENABLE && LSF_QUANTIZATION_ENABLE
    bits = g729_encode_bitstream(pkt{frameIdx});
    pkt{frameIdx}.bitstream = bits;
end 
  
  
  % perform decoder operations
  if DECODE_ENABLE
    
     [s_dec,g729.dec.state,decSignals] = g729_decode( g729.dec.state,... 
                                                      pkt{frameIdx},...
                                                      g729.consts,...
                                                      g729.dec.rom,...
                                                      g729.rom.cb); 

if DEBUG_ENABLE                                                      
     d1 = (frameIdx-1)*Nf + 1;
     d2 = d1 + Nf - 1;      
     
     % save signals
     signals.dec.u(d1:d2)      = decSignals.u;
     signals.dec.s_hat(d1:d2)  = decSignals.s_hat;
     signals.dec.s_post(d1:d2) = decSignals.s_post;
     signals.dec.out(d1:d2)    = decSignals.s_dec;
     signals.dec.residual(d1:d2)    = decSignals.residual_dec;
     signals.dec.residual_pf(d1:d2) = decSignals.residual_pf_dec;
end


  end  % DECODE enable
  
  % frameStats
  frameStats.k(:,frameIdx) = k';
  frameStats.energy(frameIdx) = norm(g729.enc.state.speechBuf); 
  frameStats.T_op(frameIdx) = T_op;  % save open loop pitch
  frameStats.wi(:,frameIdx) = wi';   % save LSF
  frameStats.r_corr(:,frameIdx) = r'; % 0th lag
end

wavwrite(yFile,8000,'encIn.wav');
wavwrite(signals.enc.s(40+(1):end),8000,'enc_s.wav');
wavwrite(signals.enc.s_hat(40+(1):end),8000,'enc_s_hat.wav');
wavwrite(signals.enc.sw(40+(1):end),8000,'enc_sw.wav');
wavwrite(signals.enc.u(40+(1):end),8000,'enc_u.wav');
wavwrite(signals.enc.ew(40+(1):end),8000,'enc_ew.wav');
wavwrite(signals.dec.s_hat(40+(1):end),8000,'dec_s_hat.wav');
wavwrite(signals.dec.s_post(40+(1):end),8000,'dec_s_post.wav');
wavwrite(signals.dec.out(40+(1):end),8000,'dec_out.wav');
