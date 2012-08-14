function [O,OFR,S,SFR,D,U,M,SE] = newfp_onsets(D,SR,B,T,blocking)
% [O,OFR,S,SFR,D,U,M,SE] = newfp_onsets(D,SR,B,T,blocking)
%   Calculate onsets for new FP system
%   D is a waveform at sampling rate SR.
%   Break into B bands, calculate onsets in each 
%   at average rate of T per second.
%   Return onset mask O (with column frame rate OFR, includes blocking)
%   and actual complex subbands S (at full frame rate SFR).
%   blocking causes envelopes to be pooled over this many
%   subsamples before onset calculation is applied.
%   D,U,M are diagnostic outputs of adaptiveonsets.
% 2010-10-10 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  B = 8;  end
if nargin < 4;  T = 4;  end
if nargin < 5;  blocking = 4;  end

deadtime = 128;

%blocking = 4; % how many subband samples to combine
overblock = 2; % how much longer than minimum is window

% Convert to 8 bands with specgram (actually 8+1 = 9)
%S = specgram(D,2*B);

use_fft = 0;

if use_fft
  % Now we use a 128 point FFT ...
  fftlen = 2*B*blocking*overblock;
  % .. and step it by 25% (32 points)
  ffthop = B*blocking;
  S = abs(specgram(D,fftlen,0,fftlen,fftlen-ffthop));
  % .. and combine energy in successive blocks of 8 FFT bins
  % to give 8 energies per hop
  for r = 1:B
    S(r,:) = sqrt(sum(S((r-1)*blocking*overblock + [1:blocking*overblock],:).^2));
  end
  S = S(1:B,:);

  SFR = SR/ffthop;

else
  % not FFT, use cosine-modulated filterbank
  assert(B==8); % because that's all Subband_Analysis_8 can do
  S = abs(Subband_Analysis_8(D')).^2;
  %S = abs(subband_analysis_helper(D)).^2;
  % Now integrate blocks of <blocking> from windows of
  % <blocking*overblock>

  hop = 4;
  nsm = 8;
  hw = hanning(nsm);
  nc = floor(size(S,2)/hop)-(floor(nsm/hop)-1);
  SS = zeros(B,nc);
  for i = 1:nc;
    SS(:,i) = sqrt(S(:,(i-1)*hop+[1:nsm])*hw);
  end
  
  S = SS;

  SFR = SR/(B*hop);

end

%size(S)

OFR = SFR;

SE = abs(S);
%size(SE)

% emphasize onsets
%hpf_pole = 0.997;
hpf_pole = 0.98;
%SE = (filter([1 -1],[1 -hpf_pole],SE')');
%SE = filter(hanning(6)',1,SE')';

bn = [0.1883 0.4230 0.3392 0 -0.3392 -0.4230 -0.1883];
%SE = filter(bn,[1 -hpf_pole],SE')';

% mean separation between onsets in steps
% Want T onsets/sec -> sep time is how many frames?
%meantsep = OFR/T;
meantsep = OFR/T;
[O,D,U,M] = adaptiveonsets(SE,meantsep,deadtime);

% to match what happens inside adaptiveonsets
SE = filter(bn,[1 -hpf_pole],SE')';
