function [R,HO] = ht_match(H,W,IX)
% [R,HO] = ht_match(H,W,IX)
%     Match a set of hashes H against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames (of 32ms).  Positive <modalDT> means the
%     query matches after the start of the reference track.
%     HO returns the actual hashes that this implies for IX'th return.
%     as rows of <time in match vid> hash <timeskew of query>
%     W is the width (in time frames) of the coarse quantization of the 
%     time differences.  Differences within this window count as matches.
% 2010-10-24 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  W = 1;  end
if nargin < 3;  IX = 1;  end

%disp([num2str(size(H,1)),' hashes']);
Rt = ht_get_hits(H);
nr = size(Rt,1);

if nr > 0

  % Find all the unique tracks referenced
  [utrks,xx] = unique(sort(Rt(:,1)),'first');
  utrkcounts = diff([xx',nr]);

  [utcvv,utcxx] = sort(utrkcounts, 'descend');
  % Keep at most 100 per hit
  utcxx = utcxx(1:min(100,length(utcxx)));
  utrkcounts = utrkcounts(utcxx);
  utrks = utrks(utcxx);
  
  nutrks = length(utrks);
  R = zeros(nutrks,4);
  
  for i = 1:nutrks
    tkR = Rt(Rt(:,1)==utrks(i),:);
    % Quantize times per window
    tkR(:,2) = round(tkR(:,2)/W);
    % Find the most popular time offset -- TO NEAREST SECOND??
    [dts,xx] = unique(sort(tkR(:,2)),'first');
    dtcounts = 1+diff([xx',size(tkR,1)]);
    [vv,xx] = max(dtcounts);
    %    [vv,xx] = sort(dtcounts, 'descend');
    %R(i,:) = [utrks(i),vv(1),dts(xx(1)),size(tkR,1)];
    % Keep everything with
    R(i,:) = [utrks(i),sum(abs(tkR(:,2)-dts(xx(1)))<=1),dts(xx(1)),size(tkR,1)];
  end

  % Sort by descending match count
  [vv,xx] = sort(R(:,2),'descend');
  R = R(xx,:);

  % Extract the actual landmarks
  % maybe just those that match time?
  %H = Rt((Rt(:,1)==R(IX,1)) & (abs(Rt(:,2)-R(IX,3))<=1),:);
  % no, return them all
  Hix = find(Rt(:,1)==R(IX,1));
  HO = [Rt(Hix,:),Rt(Hix,2)];
  % Restore the original times
  for i = 1:length(Hix)
    %HO(i,:) = [(Rt(Hix(i),:)), Rt(Hix(i),2)];
    hqix = find(H(:,2)==Rt(Hix(i),3));
    %hqix = hqix(1);  % if more than one...
    HO(i,1) = HO(i,1)+H(hqix(1),2);
  end


  % Return no more than 100 hits, and only down to 10% the #hits in
  % most popular
  maxrtns = 100;
  if size(R,1) > maxrtns
    R = R(1:maxrtns,:);
  end
  maxhits = R(1,2);
  nuffhits = R(:,2)>(0.1*maxhits);
  %R = R(nuffhits,:);

else
  R = zeros(0,4);
  disp('*** NO HITS FOUND ***');
end
