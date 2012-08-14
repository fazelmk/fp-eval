function [D, SR] = readaudio(FN,S,FM,DS)
% [D, SR] = readaudio(FN,S,FM,DS)
%   Read in an audio file, using wavread, mp3read, or m4aread as
%   appropriate. 
%   S is an optional sample range; FM optionally forces mono; 
%   DS is an optional downsample factor (by 2 or 4).
% 2010-09-16 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; S = 0; end
if nargin < 3; FM = 0; end
if nargin < 4; DS = 1; end

[pth,nam,ext] = fileparts(FN);
ext = lower(ext);
if strcmp(ext,'.mp3')
  [D,SR] = mp3read(FN,S,FM,DS);  % always mono/downsample
elseif strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4')
  [D,SR] = m4aread(FN,S,FM,DS);
else
  if S == 0
    [D,SR] = wavread(FN);
  else 
    [D,SR] = wavread(FN, S*DS);
  end
  if FM==1
    D = mean(D,2);
  end
  if DS > 1
    D = resample(D,1,DS);
    SR = SR/DS;
  end
end

