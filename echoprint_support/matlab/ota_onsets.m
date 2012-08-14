function [B,T,M] = ota_onsets(D,SR,tdens,blocking,bands)
% [B,T,M] = ota_onsets(tk) or ota_onsets(d,sr,tdens,blocking,bands)
%    Load file tk, or take waveform data (d,sr) 
%    and calculate the per-band onset times (only).
%    B and T give the band and time of each
%    onset. M gives its "prominence".
% 2011-04-17 Dan Ellis dpwe@ee.columbia.edu after calc_onset_ftrs

if nargin < 3; tdens = 1; end
if nargin < 4; blocking = 4; end
if nargin < 5; bands = 8; end

% Read in audio, and convert to mono/8 kHz
if ischar(D)
  fname = D;
  [D,SR] = readaudio(fname,0,1,2);
  if SR < 11025
    [D,SR] = readaudio(fname,0,1,1);
  end
end
%targetSR = 8000;
targetSR = 11025;
targetchans = 1;
% Fix chans/SR
[D,SR] = normaudio_ch_sr(D,SR,targetchans,targetSR);

% Whiten first?
npoles = 40;
WHITEN = 1;
if WHITEN
  D = whiten(D,npoles);
else
  disp('no whitening');
end

% Convert into onsets (in mask O at frame rate FR)
% with complex subband samples in S.
% 8 subbands of 500 Hz (maybe 9)
%bands = 8;
% 4 onsets/sec/subband (for now)
%tdens = 4;  % dropping from 4 to 2 kills d3n (17 matches -> 4)

%blocking = 4;

[On,OFR,S,SFR,DD,UU,MM] = newfp_onsets(D,SR,bands,tdens,blocking);
[nSrows, nScols] = size(S);

nrows = size(On,1);

% Pre-allocate tables to store results
tsize = 1024;
B = zeros(tsize,1);
T = zeros(tsize,1);
M = zeros(tsize,1);
specix = 0;

for b = 1:nrows
  onsettimes = find(On(b,:));
  nspecs = length(onsettimes);
  
  if specix + nspecs >= tsize
    % grow table
    B = [B;zeros(tsize,1)];
    T = [T;zeros(tsize,1)];
    M = [M;zeros(tsize,1)];
    tsize = size(S,1);
  end
  
  B(specix + [1:nspecs]) = b;
  T(specix + [1:nspecs]) = onsettimes/OFR;
  M(specix + [1:nspecs]) = MM(b,onsettimes);
  
  specix = specix + nspecs;
end

B = B(1:specix);
T = T(1:specix);
M = M(1:specix);
