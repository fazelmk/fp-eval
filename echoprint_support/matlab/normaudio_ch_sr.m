function [D,SR] = normaudio_ch_sr(D,SR,CH,TSR)
% [D,SR] = normaudio_ch_sr(D,SR,CH,TSR)
%   Normalize audio channel count and sampling rate.
%   If CH is specified and > 0, coerce to that many channels.
%   If TSR is specified and > 0, resample to that sampling rate.
% 2010-10-10 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  CH = 0;  end
if nargin < 4;  TSR = 0; end

[nr,nc] = size(D);

% Enforce requested channel count and sampling rate, if any.

if (CH > 0) && (nc ~= CH)
  if CH == 1
    D = mean(D,2);
  else
    % distribute channels one at a time, with repeats
    D = D(1+rem([0:(CH-1)],nc),:);
  end
end

if (TSR > 0) && (SR ~= TSR)
  D = resample(D,TSR,SR);
  SR = TSR;
end
