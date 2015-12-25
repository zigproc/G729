% G729 encoder ROM initialization
clear g729;

g729.consts.frameSize = 80;
g729.consts.subframeSize = 40;
g729.consts.speechBufSize = 240;
g729.consts.lpcOrder = 10;


% preproc filter
% Equation (1)
g729.enc.h1.b = [0.46363718 -0.92724706 0.46363718];
g729.enc.h1.a = [-1.9059465 0.9114024];

% LP analysis windowing
% Equation (3)
w = zeros(1,240);
w(1:200) = 0.54 - 0.46*cos(2*pi*(0:199) / 399);
w(201:end) = cos(2*pi*(0:39) / 159);

g729.enc.w_lp = w;

% autocorrelation bandwidth expansion
% Equation (6)

f0 = 60;
fs = 8000;

w = exp(-0.5*(2*pi*f0*(1:10) / fs).^2);
g729.enc.w_lag = w;

% LP to LSP conversion
% 60 equally spaced points between 0 and pi
% also include subdivided by 4 intervals
% original points are w(1:4:end)
% q_i = cos(w_i)
% w_i goes from 0 to pi, q_i goes from 1 to -1
% TBD don't need old method.  Just use bisection method

N = 60;
numPoints = 4*N - 3;
w = linspace(0,pi,numPoints);
g729.enc.lsp_qi = cos(w);

% interpolation filters for ACB search
Ninterp = 3;  % pitch search needs 1/3 resolution
fs_interp = Ninterp*fs;
fc_interp_r = 3600;
fc_interp_u = 3600;
wc_r = fc_interp_r / (fs_interp/2);   % normalized cutoff freq 
wc_u = fc_interp_u / (fs_interp/2);   

x_interp = -11:11;   % sample points for sinc
% build interpolation filter for R(k)
h_interp = [0 hamming(length(x_interp))' .* sinc(wc_r * x_interp)];
h_interp = Ninterp*(h_interp / norm(h_interp,1));
g729.enc.acb.interp_r = h_interp;  % interpolate R(k)

% build interpolation filter for u(n)
x_interp = -29:29;
h_interp = [0 hamming(length(x_interp))' .* sinc(wc_u * x_interp)];
h_interp = Ninterp*(h_interp / norm(h_interp,1));
g729.enc.acb.interp_u = h_interp;  % interpolate past excitation


% postproc filter
% Equation (TBD)
g729.dec.rom.h2.b = [0.93980581 -1.8795834 0.93980581];
g729.dec.rom.h2.a = [-1.9330735 0.93589199];

% dec needs to create Adaptive codebook vector
g729.dec.rom.interp_u = g729.enc.acb.interp_u;
g729.dec.rom.gammaN = 0.55; 
g729.dec.rom.gammaD = 0.7; 
g729.dec.rom.gammaP = 0.5; 
g729.dec.rom.gammaT(1) = 0.2;  % gammaT when k1prime is postitive
g729.dec.rom.gammaT(2) = 0.9;  % gammaT when k1prime is negative
g729.dec.rom.agcK = 0.85;  % Equation (90)

% build interpolation filter for rhat (used in long term postfilter)
x_interp = -16:16;
Ninterp = 8;  % dec postfilter needs 1/8 resolution
fc_interp_rhat = 3600;
fs_interp = Ninterp*fs;
wc_rhat = fc_interp_rhat / (fs_interp/2);   % normalized cutoff freq
h_interp = [hamming(length(x_interp))' .* sinc(wc_rhat * x_interp) zeros(1,7)]; % extend filter to 8*5 length
h_interp = Ninterp*(h_interp / norm(h_interp,1));
g729.dec.rom.interp_rhat = h_interp;

clear w x_interp h_interp wc_r wc_u wc_rhat;
clear f0 fs fs_interp Ninterp fc_interp_u fc_interp_r fc_interp_rhat N numPoints;



