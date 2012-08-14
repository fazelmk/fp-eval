function T = gentdiffs(O,B,N,C)
% T = gentdiffs(O,B,N,C)
%     O is a list of onstimes giving detected onsets; B is the
%     corresponding band index for each one.
%     T returns a set of rows 
%     of [time band tdiff ...] with N tdiffs per row, where tdiff 
%     is the time difference between successive onsets in a
%     particular band.  C controls fanout: the N tdiffs are 
%     taken as all possible choices of N succeeding onsets in 
%     the C onsets following each base onset; C>=N, and the number 
%     of onsets returned is nchoosek(C,N)*length(O) (or so).
% 2011-04-16 Dan Ellis dpwe@ee.columbia.edu

nbands = max(B);

% preallocate an output buffer (that will expand as needed)
tsize = 1024;
T = zeros(tsize, 2+N);
nts = 0;

% Generate all config patterns (subsets of each set of successive
% onset times that are used to calculate the stored time diffs).
% The result is a fixed set of indices given N and C; there's no
% need to continually recalculate it, or indeed to calculate it at
% run-time at all.
configs = nchoosekperms(C,N);

% Make first index always 1 (since we always diff against the first
% item) then make others start at 2
configs = [ones(size(configs,1),1),1+configs];

for b = 1:nbands
  % extract all the onset times for this band in ascending order
  btimes = sort(O(B==b));
  
  for t = 1:(length(btimes)-C)
    % C+1 successive onsets, i.e. base onset plus the next C
    ots = btimes(t+[0:C])';
  
    % Generate the new rows of the form
    %   bastime band timediff1 timediff2 ...
    % diff(X,1,2) takes 1st-order differences along the 2nd
    % dimension of X (i.e. returns the time differences in each row)
    newTs = [repmat([ots(1),b],size(configs,1),1), diff(ots(configs),1,2)];

    ntoadd = size(newTs,1);
    while (nts + ntoadd) > tsize
      % double the size of the array allocated to hold hashes
      T = [T;zeros(tsize,size(T,2))];
      tsize = size(T,1);
    end

    % store the new results
    T(nts + [1:ntoadd],:) = newTs;
    nts = nts + ntoadd;
  end
end

% trim the unused, pre-allocated rows
T = T(1:nts,:);
