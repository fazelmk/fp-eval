function [Y,A] = whiten(X,P,B,T)
% y = whiten(x,p,B,T)
%  Fit a p'th order LPC to the whole of X, and inverse-filter by
%  it.  B is a blocklen for online whiten (default 10000).
%  T is time const for online whiten (default 8).
% 2010-11-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  P = 40; end
if nargin < 3;  B = 10000; end
if nargin < 4;  T = 8; end

MATLAB = 0;

if MATLAB
  A = lpc(X,P);
  Y = filter(A,1,X);
else
  % not identical, but does the job.
  blocklen = B;
%  T = 8;
  Y = zeros(length(X),1);
  st = [];
  nblks = floor(length(X)/blocklen);
  for i = 1:nblks
    xx = (i-1)*blocklen + [1:blocklen];
    [Y(xx),st] = online_whiten(X(xx),st,T,P);
  end
  xx = (nblks*blocklen+1):length(X);
  if length(xx) > 0
    Y(xx) = online_whiten(X(xx),st,T,P);
  end
end


