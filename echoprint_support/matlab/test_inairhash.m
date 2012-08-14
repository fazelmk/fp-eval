% test_inairhash
%
% Test of the new dt-based inair hash
% 
% 2011-04-16 Dan Ellis dpwe@ee.columbia.edu

% Get a list of all the items

inairdir = 'inair_wednesday';

recalc_ref = 1;

%dtres = 32/11025;
%dtres = 128/11025;   % doesn't seem to help any??

% How many time differences (gaps between onsets) to encode in each
% hash.  A larger number leads to a larger hash space, but a lower 
% probability of hashes matching between true matches
%ndts = 2;
% From how many onsets following each base to choose the onsets
% whos time differences will be recorded, so must be >= ndts.  
% The number of hashes recorded grows as nchoosek(fanoutcontext, ndts)
% so you don't want it too great.  But having more hashes improves 
% the chances of some of them matching (proportionally).
%fanoutcontext_ref = 3;
%fanoutcontext_air = 3;

if recalc_ref

  fl = myls(fullfile(inairdir,'*.mp3'));

%  fl = fl(1:20);

  nfl = length(fl);

  ht_clear();

%  disp(['numdts = ', num2str(ndts), ' context_ref = ', num2str(fanoutcontext_ref)]);
%  disp([' context_air = ', num2str(fanoutcontext_air)]);

  % Put all the reference items in the database

  totnH = 0;
  totnS = 0;
  tstart = now();

  atres = 256/11025;  % units for absolute time resolution in hash
                      % for figuring the duration of each track only
  
  for i = 1:nfl
    %H = newfp_ota(fl{i},0,dtres,ndts,fanoutcontext_ref);
    H = newfp_ota(fl{i},0);
    [p,n,e] = fileparts(fl{i});
    ht_store(H,n);
    disp(['added ',num2str(length(H)),' hashes for ',n]);
    totnH = totnH + length(H);
    totnS = totnS + max(H(:,1))*atres;
  end

  atime = 3600*24*(now()-tstart);

  disp(['Added ',num2str(nfl),' tracks in ', num2str(atime), ...
       ' sec (proc time = ',num2str(atime/totnS),' s per sec of audio)']);
  disp([num2str(totnH),' hashes for ', num2str(totnS), ' s = ', ...
        num2str(totnH/totnS), ' hash/sec on avg']);

end % recalc_ref
  
% build recognition matrix
reccounts = zeros(nfl, nfl);

tstart = now();

qts = 0;

for i = 1:nfl;
  % hashes for inair version
  [p,n,e] = fileparts(fl{i});
  %Hn = newfp_ota(fullfile(p,[n,'.wav']),0,dtres,ndts,fanoutcontext_air);
  Hn = newfp_ota(fullfile(p,[n,'.wav']));
  fn{i} = n;
  R = ht_match(Hn);
  reccounts(i,R(:,1)) = R(:,2);
  %qts = qts + length(dn)/srn;
  qts = qts + max(Hn(:,1))*atres;
end

subplot(121)
imagesc(log10(1+reccounts));
colorbar;
%title(['ndts=',num2str(ndts),' ctxt_air=',num2str(fanoutcontext_air), ...
%       ' - log10(1+reccounts)']);
set(gca,'YTick',[1:nfl]);
set(gca,'YTickLabel',fn);

subplot(222)
tru = diag(reccounts);
fls = max(reccounts'-diag(diag(reccounts')))';
nit = length(tru);
plot(1:nit,tru,1:nit,fls,'-r');
grid
axis([0 nit+1 0 100])
title('blue = true count; red = worst false count');

subplot(224)
plot(1:nit, tru./fls);
grid
axis([0 nit+1 0 10])
title('true count / worst false count')

qtime = 3600*24*(now()-tstart);

disp(['Ran ',num2str(nfl),' queries in ',num2str(qtime), ...
      ' sec (proc time = ',num2str(qtime/qts),' s / sec audio)']);
tth = 20;
disp([num2str(sum(tru<tth)),' / ',num2str(nfl), ...
      ' true matches with < ', num2str(tth),' hits']);
rth = 2;
disp([num2str(sum((tru./fls)<rth)),' / ',num2str(nfl), ...
      ' true/false ratios < ', num2str(rth)]);
