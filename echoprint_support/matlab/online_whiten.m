function [y, s_o] = online_whiten(x, s_i, t_dec, order)
% [y, s_o] = online_whiten(x, s_i, alpha, order)
% - get a new block of samples
% - update autocorrelation for this new block
% - newly calculate inverse filter
% - apply inverse filter
%   - crossfade with previous coefficients?
% - return new block of filtered samples

if nargin < 2;  s_i = []; end
if nargin < 3;  t_dec = 100; end
if nargin < 4;  order = 8; end

mexfname = 'online_whiten_helper';
HAVE_MEX = (exist(mexfname)==3);
if HAVE_MEX == 0
  try 
    disp(['compiling ',mexfname,'...']);
    mex([mexfname,'.c']);
  end
end
HAVE_MEX = (exist(mexfname)==3);

%HAVE_MEX = 0;

if HAVE_MEX

  % reset on first call
  if length(s_i) == 0
    online_whiten_helper([]);
  end
  
  y = online_whiten_helper(x, order, t_dec);

  % Mark that future calls are not first call
  s_o = 1;

else

  % Matlab version
  
  if length(s_i) == 0
    state.order = order;
    state.autoco = zeros(order+1, 1);
    state.hist_in = zeros(order, 1);
    state.alpha = 1/t_dec;
%  state.lasta_i = zeros(1, order+1);
%  state.lasta_i(1) = 1;
  else
    state = s_i;
  end

  N = length(x);

  % calculate autocorrelation of current block
  aco = zeros(order+1, 1);
  for i = 0:order
    for j = 0:(N-1)
      if (j-i) < 0
        %aco(1+ i) = aco(1+ i) + x(1+ j) * state.hist_in(1+ order + j-i);
      else
        aco(1+ i) = aco(1+ i) + x(1+ j) * x(1+ j-i);
      end
    end
    % normalize
    %aco(1+ i) = aco(1+ i)/N;
  end

  % update autocorrelation
  state.autoco = state.autoco + state.alpha*(aco - state.autoco);

  % calculate new filter coefficients
  a_i = zeros(1,order+1);
  a_i(1)  = 1;

  % Durbin's recursion, per p. 411 of Rabiner & Schafer 1978
  R = state.autoco;

  E = R(1+ 0);
  p = order;
  %alphaM1 = zeros(1,p);
  alpha = zeros(1,p);
  for i = 1:p
    sumalphaR = 0;
    for j = 1:i-1
  %    sumalphaR = sumalphaR + alphaM1(j)*R(1+ i-j);
      sumalphaR = sumalphaR + alpha(j)*R(1+ i-j);
    end
    ki = (R(1+ i) - sumalphaR)/E;
    alpha(i) = ki;
  %  for j = 1:i-1
  %    alpha(j) = alphaM1(j) - k(i)*alphaM1(i-j);
  %  end
    for j = 1:floor(i/2)
      aj = alpha(j);
      aimj = alpha(i-j);
      alpha(j) = aj - ki*aimj;
      alpha(i-j) = aimj - ki*aj;
    end
    E = (1-ki^2)*E;
  %  alphaM1 = alpha;
  end
  a_i(2:(p+1)) = -alpha;

  %R
  %a_i

  % caculate new output
  y = zeros(N,1);

  for i = 0:N-1
    acc1 = x(1+ i);
  %  acc2 = x(1+ i);
    for j = 1:order
      if i-j < 0
        x_o = state.hist_in(1+ order + i-j);
        acc1 = acc1 + a_i(1+ j)*x_o;
  %      acc2 = acc2 + state.lasta_i(1+ j)*x(1+ i-j);
      else
        acc1 = acc1 + a_i(1+ j)*x(1+ i-j);
  %      acc2 = acc2 + state.lasta_i(1+ j)*x(1+ i-j);
      end
    end
  %  y(1+ i) = (i/N)*acc1 + ((N-i)/N)*acc2;
    y(1+ i) = acc1;
  end

  % save state
  state.hist_in = x(N-[p:-1:1]);
  %state.lasta_i = a_i;

  s_o = state;

end
