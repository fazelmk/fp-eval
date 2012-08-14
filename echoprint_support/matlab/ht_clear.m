function ht_clear()
% ht_clear()
%  Access the persistent store to reset the hash table.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

global HashTable HashTableCounts HashTableNames HT_Hsize HT_Rsize

% Make sure the hash mex is set
mexfname = 'jenkinshash';
HAVE_MEX = (exist(mexfname)==3);
if HAVE_MEX == 0
  try 
    disp(['compiling ',mexfname,'...']);
    mex([mexfname,'.c']);
  end
end
HAVE_MEX = (exist(mexfname)==3);


%if exist('HashTable','var') == 0
%   HashTable = [];
%end

nhashes = 2^20;

% 1M hashes x 32 bit entries x 100 entries = 400MB in core
%maxnentries = 100;
maxnentries = 100;

disp(['Max entries per hash = ',num2str(maxnentries)]);

%if length(HashTable) == 0
  HashTable = zeros(maxnentries, nhashes, 'uint32');
  HashTableCounts = zeros(1, nhashes);
%end

% Reset the table that maps hash table indices to names
HashTableNames = [];

% Reset stats on numbers of hashes queried and hits returned
HT_Hsize = 0;
HT_Rsize = 0;

