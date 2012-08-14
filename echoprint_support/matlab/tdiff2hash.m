function H = tdiff2hash(L, dtres, atres)
% H = tdiff2hash(L, dtres, atres)
%    Convert time-difference landmarks to hashes.
%    Rows of L are [time band dt1secs [dt2secs ...]]
%    dtres is the time quantization step size for time differences (32/11025)
%    atres is the time quantization step size for absolute times (256/11025)
% 2011-04-16 Dan Ellis dpwe@ee.columbia.edu

% Note: newfp_ota passes in dtres (but not atres), so the default isn't used 
if nargin < 2;  dtres = 128/11025; end  % default 11.6ms time diff quantization
if nargin < 3;  atres = 256/11025; end  % default 23.2ms abs time quantization

%tres = newfp_time_res();

T = round(L(:,1)/atres);
B = L(:,2);

dts = round(L(:,3:end)/dtres);

% allocate 10 bits for each time diff
hsz = 2^10;

maxb = 9;

ndtcols = size(dts,2);
dtqval = dts * (hsz.^[0:ndtcols-1]');

H = [T, B + maxb*dtqval];

