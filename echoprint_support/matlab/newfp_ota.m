function H = newfp_ota(D,SR,dtres,ndts,fanoutcontext)
% H = newfp_ota(D,SR,dtres,ndts,fanoutcontext)
%  Calculate the new over-the-air capabale subband event hashes.
%  D,SR are a sound waveform, or D is a sound file name.
%  Load, normalize, process, return a sequence of 
%  <time hash> fingerprints.
%  Optional ndts is the number of time differences to use in each
%  hash (default 2); fanoutcontext is the number of succeeding
%  events to choose these differences from (default 4).
% 2011-04-17 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; dtres = 256/11025; end
% How many time differences to include
if nargin < 4; ndts = 2; end
% .. out of how many succeeding events?
if nargin < 5; fanoutcontext = 4; end

% Maybe read in a file if we're passed a name
if ischar(D)
  fname = D;
  [D,SR] = readaudio(fname);
end

% Calculate all the subband onset times
[bn,tn,mn] = ota_onsets(D,SR);
% Generate the fanned-out sets of time differences, and hash them
H = tdiff2hash(gentdiffs(tn,bn,ndts,fanoutcontext), dtres);
