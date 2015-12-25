function [s_post,rhat_curr,rhat_pf,decState] = dec_postfilter( decState,a_hat,Nsf,T,gammaN,gammaD,gammaP,interp_r)
  
  ENABLE_TILT_COMPENSATION = 1;
  ENABLE_AGC = 1;
  
  M = length(a_hat);  
  assert(length(decState.rhat) >= decState.ZI_RHAT + Nsf -1);
  assert(length(decState.s_hat) >= decState.ZI_S_HAT + Nsf - 1);
  assert(decState.ZI_RHAT - M >= 1);
  assert(decState.ZI_S_HAT - M >= 1);
  assert(decState.ZI_RHAT - (T+4) >= 1);
   
  ZI_RHAT = decState.ZI_RHAT;
  ZI_S_HAT = decState.ZI_S_HAT;
    
  % sfPrime in spec
  s_post = zeros(1,Nsf);
     
  % initalize state for computing rhat, which is past input signal s_hat
  fb_state = decState.s_hat(ZI_S_HAT - M: ZI_S_HAT - 1);
  
  % compute residual (hat) by filtering reconstructed speech 
  % thru prediction error filter A_hat(z/gammaN) = numerator of Short term postfilter
  % Equation (79),(84)
  
  [decState.rhat(ZI_RHAT:ZI_RHAT+Nsf-1),~] = allzero_filter(decState.s_hat(ZI_S_HAT:ZI_S_HAT+Nsf-1),...
                                                            fb_state,...
                                                            a_hat,...
                                                            gammaN);
  
  rhat    = decState.rhat;
  rhat_curr = rhat(ZI_RHAT:ZI_RHAT+Nsf-1);
  
  rhat_pf = long_term_postfilter(rhat,...
                                 ZI_RHAT,...
                                 interp_r,...
                                 T,...
                                 gammaP,...
                                 Nsf);
  
  % filter rhat thru short term synthesis filter 1/(gf*Ahat(z/gammaD))
  % -------------------------------------------------------------------
  [gf,hf] = compute_stp_gain(a_hat,...
                             gammaN,...
                             gammaD);
                        
  % synthesis filtering of rhat_pf
  [s_hat_pf_syn,decState.stp.ahat.fb ] = allpole_filter(rhat_pf,...
                                                        decState.stp.ahat.fb ,...
                                                        a_hat,...
                                                        gammaD);
  s_hat_pf_syn = (1/gf)*s_hat_pf_syn;
  
  if( ENABLE_TILT_COMPENSATION)
    % filter output of synthesis filter thru tilt compensation filter to get sf(n)
    [sf,decState.htilt.ff] = tilt_compensation(s_hat_pf_syn,decState.htilt.ff,hf);
  else  
    sf = s_hat_pf_syn;
  end
  
  
  if( ENABLE_AGC)
    % scale sf(n) using adaptive gain control
    s1 = decState.s_hat(ZI_S_HAT:ZI_S_HAT + Nsf - 1);
    s2 = sf(1:Nsf);
    Es2 = s2*s2';
    
    G = 0;
    if( Es2 ~= 0)
      G = (s1*s1') / (s2*s2');
    end
    
    % Equation (90)
    for n=0:Nsf-1
      gn = 0.85*decState.gn_prev + 0.15*G;
      s_post(n+(1)) = gn*sf(n+(1));
      decState.gn_prev = gn;
    end
    
  else
    s_post = sf;
  end
    
end


function [y,state] = allzero_filter(s,state,a_hat,gamma)

  assert(length(state) == length(a_hat));
  
  y = zeros(1,length(s));
  
  a_hat = fliplr(a_hat .* (gamma .^ (1:length(a_hat))));
  
  for n=0:(length(s)-1)
    sample = s(n+(1));
    sum    = sample;
    sum   += (a_hat*state'); 
    
    state = [state(2:end) sample];
    
    y(n+(1)) = sum;
  end  
end


function [y,state] = allpole_filter(s,state,a_hat,gamma)

  assert(length(state) == length(a_hat));
  
  y = zeros(1,length(s));
  
  a_hat = fliplr(a_hat .* (gamma .^ (1:length(a_hat))));
  
  for n=0:(length(s)-1)
    sample = s(n+(1));
    sum  = sample;
    sum -= a_hat*state';
    
    state = [state(2:end) sum];
    y(n+(1)) = sum;
  
  end
end

function [y,state] = tilt_compensation(s,state,h)
  
  assert(length(h) == 20);
  assert(length(s) == 40);
  
  y = zeros(1,40);
  % Equation (87)
  R0 = h*h';
  R1 = h(1:19)*h(2:20)';
  k1Prime = -R1/R0;
  
  if(k1Prime >=0)
    gammaT = 0.2;
  else  
    gammaT = 0.9;
  end
  
  G = gammaT*k1Prime;
  
  gt = 1 - abs(G);
  
  % filter S using Ht(z) = (1/gt)*(1+Gz^-1)
  % Equation (86)
  
  for n=0:39
    sample = s(n+(1));
    sum = sample + G*state;
    state = sample;
    y(n+(1)) = (1/gt)*sum;
  end
end

function [gf,hf] = compute_stp_gain(a_hat,gammaN,gammaD)

  % Equation (85)
  assert(length(a_hat) == 10);
  
  a_hat_N = a_hat .* (gammaN .^ (1:10));   % treat as input signal
  a_hat_D = fliplr(a_hat .* (gammaD .^ (1:10)));  % treat as filter
  hf = zeros(1,20);
  
  % compute truncated impuse response of Hf(z), n=0..19
  in = [1 a_hat_N zeros(1,9)];
  state = zeros(1,10);
  
  % run input thru denominator of STP, given by Ahat(z/gammaD)
  hf = allpole_filter(in,state,a_hat,gammaD);
  
  % compute 1-norm
  gf = norm(hf,1);  

end


function rhat_pf = long_term_postfilter(rhat,ZI_RHAT,interp_r,T,gammaP,Nsf)
  
  R       = zeros(1,3);
  rn      = rhat(ZI_RHAT:ZI_RHAT+Nsf-1);
  Erhat   = rn*rn';
  
  % compute correlation in range [T-1,T+1]
  % Equation (80)
  kMin = (T-1);
  kMax = (T+1);
  for k=kMin:kMax
    rk = rhat(ZI_RHAT-k:ZI_RHAT-k+Nsf-1);  %delayed signal at integer resolution
    R(k-kMin+(1)) = rn*rk';
  end
  
  [mx,ix] = max(R);
  Topt = kMin + ix - (1);
  
  % phases of interpolation filter
  assert(mod(length(interp_r),4)==0);
  bPhase{1} = fliplr(interp_r(1:8:end));
  bPhase{2} = fliplr(interp_r(2:8:end));
  bPhase{3} = fliplr(interp_r(3:8:end));
  bPhase{4} = fliplr(interp_r(4:8:end));
  bPhase{5} = fliplr(interp_r(5:8:end));
  bPhase{6} = fliplr(interp_r(6:8:end));
  bPhase{7} = fliplr(interp_r(7:8:end));
  bPhase{8} = fliplr(interp_r(8:8:end));
 
  % interpolate in 1/8 resolution around Topt.  Not interpolating
  % R values, but interpolating rhat(n-k) at various fractional resolution k 
  % to compute the correlation
  
  % start at -Topt-1, and scan with 1/8 resolution up to -Topt+1
  % total of 17 options (e.g k = [-58,-57 7/8, -57 6/8, ... -57, -56 7/8 ... -56]
  
  Rfrac = zeros(1,17);
  glk   = zeros(1,17);
  
  ks = -8*(Topt+1);
  for i=0:16
    kInt = floor((ks+i)/8);
    kPhase = mod(ks+i,8);
  
    [Rfrac(i+(1)),glk(i+(1))] = compute_R(rhat,ZI_RHAT,-kInt,bPhase{kPhase+(1)},Nsf);
  end
                                            
  [maxCorr,ix] = max(Rfrac);
  
  % check prediction gain
  % Equation (82)
  if( Erhat ~= 0 && (maxCorr^2 / Erhat) < 0.5)
    gl = 0;
  else
    gl = glk(ix);
    gl = max(gl,0.0);
    gl = min(gl,1.0);
  end
  
  % compute rnk for optimal fractional delayed 
  ksOpt  = ks + ix - (1);
  kInt   = floor(ksOpt / 8);
  kPhase = mod(ksOpt,8);
  
  % optimal fractionally delayed signal.  Redundant, since we computed this signal
  % at some point in the search above
  rnk_opt = compute_rnk(rhat,ZI_RHAT,-kInt,bPhase{kPhase+(1)},Nsf);  
  
  % compute rhat signal after running thru long term postfilter Hp(z)
  %   rhat_pf(n) = K1*rhat(n) + K2*rhat(n-kopt) 
  G  = gammaP*gl;
  K1 = 1 / (1+G);
  K2 = G / (1+G);

  % Equation (78)
  rhat_pf = K1*rn + K2*rnk_opt;  
end

function [y,glk] = compute_R(rhat,ZI_RHAT,k,b,Nsf)


  rn  = rhat(ZI_RHAT:ZI_RHAT+Nsf-1);  % r(n) signal, undelayed
  rnk = compute_rnk(rhat,ZI_RHAT,k,b,Nsf);
  y   = 0;
  glk = 0;
  
  Enk = rnk*rnk';
  % normaized correlation 
  % Equation (81)
  
  if( Enk ~= 0)
    y = rn*rnk' / sqrt(Enk);
    glk = rn*rnk' / Enk;
  end  
end

function rnk = compute_rnk(rhat,ZI_RHAT,k,b,Nsf)

  Nw = floor(length(b)/2); % half width of filter

  rnk = zeros(1,Nsf);  % r(n-k)
  % compute r_hat(n-k) at fractional resolution
  for n=0:Nsf-1
    % conv with interpolation filter centered around rhat(n-k)
    rnk(n+(1)) = rhat(ZI_RHAT+n-k-Nw : ZI_RHAT+n-k+Nw) * b';
  end

end
