function N = ht_store(H,A)
% N = ht_store(H,A)
%   Record the set of hashes that are rows of H in persistent
%   database, associated with filename A.
%   Format of H rows are 2 columns:
%   <start_time_index> <hash>
%   time_index is ? 12 bit (i.e. integer out to 16384)
% Hash is a 32 bit int val, is hashed and masked to 20 bits
% N returns the actual number of hashes saved (excluding table overflows).
%
% 2008-12-24 Dan Ellis dpwe@ee.columbia.edu

% This version uses an in-memory global with one row per hash
% value, and a series of song ID / time ID entries per hash

global HashTable HashTableCounts HashTableNames

%if exist('HashTable','var') == 0 || length(HashTable) == 0
%   ht_clear
%end

% Fill in next slot in HashTable record
ID = length(HashTableNames) + 1;
HashTableNames{ID} = A;

maxnentries = size(HashTable,1);

nhash = size(H,1);

% Mask the hash val at 20 bits, hash it up first
Hvals = bitand(jenkinshash(uint64(H(:,2))), uint32((2^20)-1));

N = 0;

TIMESIZE = 16384;

for i=1:nhash
  toffs = mod(round(H(i,1)), TIMESIZE);
  hash = 1+Hvals(i); % avoid hash == 0
  htcol = HashTable(:,hash);
  nentries =  HashTableCounts(hash) + 1;
  if nentries <= maxnentries
	% put entry in next available slot
	r = nentries;
  else
    % choose a slot at random; will only be stored if it falls into
    % the first maxnentries slots (whereupon it will replace an older 
    % value).  This approach guarantees that all values we try to store
    % under this hash will have an equal chance of being retained.
    r = ceil(nentries*rand(1));
  end
  if r <= maxnentries
    hashval = int32(ID*TIMESIZE + toffs);
%    disp(num2str(floor(double(hashval)/TIMESIZE)));
    HashTable(r,hash) = hashval;
    N = N+1;
  end
HashTableCounts(hash) = nentries;
end
