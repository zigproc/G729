% encoder state

% speechBuf contains 15ms "previous speech", 10 ms "current frame", 5 ms lookahead
g729.enc.state.speechBuf = zeros(1,240);
% n=0 index of speech signal is 121, since keep 15 ms "prev frame" samples
g729.enc.state.ZI_S = 120+1;  

% weighted speech buf contains 143 past samples plus 80 current samples
% weighted speech history
% sw[-143] to sw[79]
g729.enc.state.weightedSpeechBuf = zeros(1,143+80);

% n=0 index of weighted speech signal is 144, since keep 143 prev samples
g729.enc.state.ZI_WS = 143+1;

% excitation buf contains 143 previous samples, plus 40 current samples, plus 9 for interpolation
% at max delay
g729.enc.state.u = zeros(1,143+40+9);
g729.enc.state.ZI_U = 143+1+9;

% initial state of prev LSP
g729.enc.state.qi_prev = cos((1:10)*pi/11);

% initial state of prev quantized LSP
g729.enc.state.qi_hat_prev = cos((1:10)*pi/11);

% initial state of previous quanitized LSF vector (used for MA prediction)
g729.enc.state.li_hat_prev = repmat((1:10)'*pi/11,1,4);

% initial state of prev "flat" state
g729.enc.state.flat_prev = 0;

% initial state of LAR
g729.enc.state.LAR_prev = zeros(1,2);

% initial state of prev quantized ACB gain
g729.enc.state.gp_hat_prev = 0.8;

% initial state of Uhat (previous quantized gain prediction errors = 20*log10(gamma)
g729.enc.state.U_hat_prev = -14*ones(1,4);

% initial states for weighted synthesis filter
% used in computation of target signal
% updated when excitation is known
g729.enc.state.wsf.w.ff = zeros(1,10); % FF state of W(z) = A(z/gamma1)/A(z/gamma2), temp random signal
g729.enc.state.wsf.w.fb = zeros(1,10); % FB state of W(z)
g729.enc.state.wsf.ahat.fb = zeros(1,10); % FB state of 1/A_hat(z)

% states used for preproc filter
g729.enc.state.preproc.h1.ff = zeros(1,2);  % FF state of H_h1(z)
g729.enc.state.preproc.h1.fb = zeros(1,2);  % FB state of H_h1(z)

%initial states for synthesis filter (locally computing synthetic speech)
g729.enc.state.sf.ahat.fb = zeros(1,10);  % FB state of 1/A_hat(z)

