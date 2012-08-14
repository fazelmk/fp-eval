function [O,D,U,M] = adaptiveonsets(E,T,deadtime)
% [O,D,U,M] = adaptiveonsets(E,T,deadtime)
%    E is an envelope (or matrix with each row an envelope).  
%    Use a diode-pump-type envelope follower 
%    to find onsets.  Adapt decay rate of diode-pump to achieve 
%    approximately 1/T onsets per sample (i.e. an average of T 
%    steps between onsets).
%    For debugging, D returns actual diode-pump level at every
%    step.  U returns the decay time at every step.
%    deadtime sets the minimum gap between onsets; this many steps
%    prior to each reported onset are cleared to zero.
%    M returns the prominence of each element in O.
% 2010-10-02

if nargin < 3
  deadtime = 32;
end

%disp(['adaptiveonsets: T=',num2str(T)])

[nr,nc] = size(E);

O = zeros(nr, nc);

mexfname = 'adaptiveonsets_helper';
HAVE_MEX = (exist(mexfname)==3);
if HAVE_MEX == 0
  try 
    disp(['compiling ',mexfname,'...']);
    mex([mexfname,'.c']);
  end
end
HAVE_MEX = (exist(mexfname)==3);

%HAVE_MEX = 1;

if HAVE_MEX
  [O,D,U,M] = adaptiveonsets_helper(E,T,deadtime);
else
  disp('Unable to use MEX optimized adaptiveonsets_helper');
  
  % Runs sequentially

  % Current decay time
  tau = ones(nr, 1);
  % time since last onset in each channel
  tsince = zeros(nr, 1);

  % Debug traces
  D = zeros(nr,nc);  % the diode pump level
  U = zeros(nr,nc);  % the decay time
  M = zeros(nr,nc);  % the "prominence" of each onset

  % Amount by which threshold exceeds actual level
  myeps1 = 0.05;

  % Current threshold
  H = zeros(nr,1); %(1+myeps1)*E(:,1);
  % ensure "contact" in first step

  % Last prominence record
  N = zeros(nr,1);
  % Previous contact status
  lcontact = zeros(nr, 1);
  
  for i = 1:nc
    contact = (E(:,i)>H);
    abv = find(contact);
    blw = find(~contact);
    attach = find( contact & (~lcontact) );
    % Remember where the threshold was when we first touched
    % but if N is not zero, keep the older (presumably smaller).
    % (N is reset to zero after deadtime)
    N(attach) = ((N(attach)==0).*H(attach))+((N(attach)>0).*N(attach));
    % Update mask
    H(abv) = E(abv,i)*(1+myeps1);
    % Apply decays
    H(blw) = H(blw) .* exp(-1./tau(blw));
    % Decays when evelope detaches
    detach = find( (~contact) & lcontact );
    if length(detach)>0
      O(detach,i) = 1;
      M(detach,i) = min(20,E(detach,i-1)./N(detach));
      % Clear any recent onset flags
      O(detach,[max(1,i-deadtime):(i-1)]) = 0;
      M(detach,[max(1,i-deadtime):(i-1)]) = 0;
      % Reset detach time counter
      tsince(detach) = 0;
    end
    % increment time-since-last-detach counters
    tsince = tsince + 1;
    % update adaptive decay thresh - increment if last onset too
    % recent, decrement if too distant
    tau = max(1, tau + ((2*(tsince<T))-1));
    % clear onset level record if it's beyond deadtime
    N(find((~contact) & (tsince>deadtime))) = 0;

    % record current contact map
    lcontact = contact;

    % debug
    D(:,i) = N; %H;
    U(:,i) = contact;
  end

end
