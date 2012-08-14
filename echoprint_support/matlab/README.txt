README for echonewfp-2011-06-04

This is the sixth code release for the new Echo Nest fingerprinter.  
It uses a new method for subband analysis, modeled on the modulated
cosine filterbank of MPEG-1 Audio.  These modified subband envelopes 
afford good recognition of over-the-air recordings, while preserving 
high accuracy for clean waveforms, at about 25 hashes/sec.

As before, the basic usage is:

>> F = newfp_ota(filename);

.. which will load an audio file and convert it to a two-column 
set of fingerprints of rows
<time hash>
where time is in seconds and hash is a positive integer smaller than 2^32.  

Putting the tracks into an actual hash table works as:

>> tks = myls('refitems/*.mp3');
>> ht_clear
>> for i = 1:length(tks); 
     ht_store(newfp_ota(tks{i}), tks{i}); 
   end

Then building a matrix of match counts is something like:

>> qs = myls('qitems/*.wav');
>> counts = zeros(length(qs), length(tks));
>> for i = 1:length(qs); 
     R = ht_match(newfp_ota(qs{i})); 
     counts(i,R(:,1)) = R(:,2);
   end


This release uses an updated MEX implemetations of adaptiveonsets 
in adaptiveonsets_helper.c, and the bit mask hash calculation 
in jenkinshash.c, and for the new subband analysis in 
subband_analysis_helper.c.

Also included is test_inairhash.m, a script which runs the entire 
system on a directory containing reference tracks as mp3 files, 
and corresponding queries as equivalently-named wav files.

The routines used are:
newfp_ota calls
  ota_onsets  which finds all the subband onsets by calling
    readaudio        to read in an audio file using mp3read
    normaudio_ch_sr  to convert it to 11 kHz/mono if not already
    newfp_onsets     to find all the onsets in the subbands which calls
      Subband_Analysis_8  to calculate the 8-subband analysis (includes compilaton of subband_analysis_helper.c)
      adaptiveonsets to return about 4 onsets per second in each band (includes compilation of adaptiveonsets_helper.c)
  gentdiffs          to map the per-band onset times into fanned-out triples
  tdiff2hash         to convert the <band tdiff1 tdiff2> into a quantized hash

Other routines are:
  ht_clear.m, ht_store.m, ht_get_hits.m, ht_match.m
                             Implementation of simple in-memory hash table.
			     Uses jenkinshash.c MEX file to hash ints.

 * end *
