

%%%%%%%%%%%%%%%  
%  DEC STATE
%%%%%%%%%%%%%%%

% initial state of prev quantized LSP
g729.dec.state.qi_hat_prev = cos((1:10)*pi/11);

% states used for postproc 2nd order filter
g729.dec.state.postproc.h2.ff = zeros(1,2);  % FF state of H_h2(z)
g729.dec.state.postproc.h2.fb = zeros(1,2);  % FB state of H_h2(z)

%initial states for synthesis filter ( computing synthetic speech)
g729.dec.state.sf.ahat.fb = zeros(1,10);  % FB state of 1/A_hat(z)

% initial states for short term postfilter denom (synthesis filtering of rhat)
g729.dec.state.stp.ahat.fb = zeros(1,10);  % FB state of A_hat(z/gammmD)
% excitation buf contains 143 previous samples, plus 40 current samples, plus 9 for interpolation
% at max delay
g729.dec.state.u = zeros(1,143+40+9);
g729.dec.state.ZI_U = 143+1+9;

% residual buf contains 143 previous samples, plus 40 current samples, plus 4 extra for interpolation
g729.dec.state.rhat = zeros(1,143+4+40);
g729.dec.state.ZI_RHAT = 143+1+4;

% reconstructed speech contains 40 current samples plus 10 previous samples
g729.dec.state.s_hat = zeros(1,10+40);
g729.dec.state.ZI_S_HAT = 10+1;

% initial state of previous quanitized LSF vector (used for MA prediction)
g729.dec.state.li_hat_prev = repmat((1:10)'*pi/11,1,4);

% initial state of Uhat (previous quantized gain prediction errors = 20*log10(gamma)
g729.dec.state.U_hat_prev = -14*ones(1,4);

% initial state of prev quantized ACB gain
g729.dec.state.gp_hat_prev = 0.8;

% initial state of prev adpative gain control (postproc), Equation (90)
g729.dec.state.gn_prev = 1.0;

% initail state for tilt compensation filter (single tap)
g729.dec.state.htilt.ff = 0;