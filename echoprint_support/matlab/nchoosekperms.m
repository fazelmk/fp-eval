function C = nchoosekperms(N,K)
% C = nchoosekperms(N,K)
%   Return all possible subsets consisting of K items from a list
%   of N items where ordering is not important.  Each row of C is a
%   different subset of 1:N (or the items in N if a list), sorted
%   into ascending order.
% 2011-04-17 Dan Ellis dpwe@ee.columbia.edu

if length(N) > 1
  L = N;
  N = length(L);
else
  L = 1:N;
end

%C = zeros(nchoosek(N,K),K);
% all possible orderings of the next N items
if N > 1
  pp = perms(L);
else
  pp = L
end

% take just the first K of them (pp(:,1:K)) 
% (i.e. sets of K items from N)
% sort them into ascending order to make permutations equivalent
% then remove the duplicate rows with unique
C = unique(sort(pp(:,1:K),2),'rows');
% we end up with nchoosek(N,K) rows of K indices in range 1:N, 
% always in ascending order, all unique.
