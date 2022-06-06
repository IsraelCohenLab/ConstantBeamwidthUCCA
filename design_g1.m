function Hd = design_g1(F1,A1,Fs,N)
%DESIGN_G1 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 05-Jan-2012 13:03:37
%

%N  = 64;     % Order
B  = 1;      %
% Fs = 16000;  % Sampling Frequency
% 
% % First Frequency Vector
% F1 = [0 0.0344827586206897 0.0689655172413793 0.103448275862069 ...
%       0.137931034482759 0.172413793103448 0.206896551724138 ...
%       0.241379310344828 0.275862068965517 0.310344827586207 ...
%       0.344827586206897 0.379310344827586 0.413793103448276 ...
%       0.448275862068966 0.482758620689655 0.517241379310345 ...
%       0.551724137931034 0.586206896551724 0.620689655172414 ...
%       0.655172413793103 0.689655172413793 0.724137931034483 ...
%       0.758620689655172 0.793103448275862 0.827586206896552 ...
%       0.862068965517241 0.896551724137931 0.931034482758621 ...
%       0.96551724137931 1];
% 
% % First Amplitude Vector
% A1 = [1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0];

h = fdesign.arbmag('N,B,F,A',N, B, F1, A1, Fs);
%h = fdesign.arbmag('N,B,F,A',4* N, B, F1, A1, Fs);
%h = fdesign.arbmagnphase('N,F,H', N, F1, A1, Fs);

Hd = design(h, 'firls');

%Hd = design(h, 'equiripple','MinPhase',true);
%Hd = design(h, 'equiripple');


% WEIGHTS=ones(N,1);
% WEIGHTS(16:end)=0.1;
% Hd = design(h, 'firls','B1Weights', WEIGHTS);


% set(Hd, 'Arithmetic', 'fixed', ...
%     'InputWordLength', 16, ...
%     'InputFracLength', 15, ...
%     'CoeffWordLength', 16, ...
%     'CoeffAutoScale', true, ...
%     'FilterInternals', 'Fullprecision');

 set(Hd, 'Arithmetic', 'single');
% set(Hd, 'Arithmetic', 'double');


% [EOF]