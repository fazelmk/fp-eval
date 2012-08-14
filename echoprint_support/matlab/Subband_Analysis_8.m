function U = Subband_Analysis_8(X)
% U = Subband_Analysis_8(X)
%    Calculate 8-band subband analysis
%    based on the MPEG-Audio Layer 1 filterbank
%    X is the input sound waveform
%    U is returned as 8 rows of subband, decimated samples.
%
% 2011-06-02 Dan Ellis dpwe@ee.columbia.edu
% based on Matlab_MPEG Subband_Analysis

mexfname = 'subband_analysis_helper';
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
  U = subband_analysis_helper(X);
else
  disp(['Unable to use MEX optimized ',mexfname]);

  % 128pt, 1/8th band low-pass prototype subsampled from Table_analysis_window
  C=[4.7700000e-07   9.5400000e-07   1.4310000e-06   2.3840000e-06   3.8150000e-06   6.1990000e-06   9.0600000e-06   1.3828000e-05   1.9550000e-05   2.7657000e-05   3.7670000e-05   4.9591000e-05   6.2943000e-05   7.6771000e-05   9.0599000e-05   1.0156600e-04  -1.0824200e-04  -1.0681200e-04  -9.5367000e-05  -6.9618000e-05  -2.7180000e-05   3.4332000e-05   1.1634800e-04   2.1886800e-04   3.3903100e-04   4.7254600e-04   6.1178200e-04   7.4720400e-04   8.6641300e-04   9.5415100e-04   9.9420500e-04   9.7131700e-04  -8.6879700e-04  -6.7424800e-04  -3.7860900e-04   2.1458000e-05   5.2213700e-04   1.1110310e-03   1.7666820e-03   2.4571420e-03   3.1418800e-03   3.7717820e-03   4.2905810e-03   4.6381950e-03   4.7521590e-03   4.5738220e-03   4.0493010e-03  3.1347270e-03  -1.8005370e-03  -3.3379000e-05   2.1615030e-03   4.7564510e-03   7.7033040e-03   1.0933399e-02   1.4358521e-02   1.7876148e-02   2.1372318e-02   2.4725437e-02   2.7815342e-02   3.0526638e-02   3.2754898e-02   3.4412861e-02   3.5435200e-02   3.5780907e-02  -3.5435200e-02  -3.4412861e-02  -3.2754898e-02  -3.0526638e-02  -2.7815342e-02  -2.4725437e-02  -2.1372318e-02  -1.7876148e-02  -1.4358521e-02  -1.0933399e-02  -7.7033040e-03  -4.7564510e-03  -2.1615030e-03   3.3379000e-05   1.8005370e-03   3.1347270e-03  -4.0493010e-03  -4.5738220e-03  -4.7521590e-03  -4.6381950e-03  -4.2905810e-03  -3.7717820e-03  -3.1418800e-03  -2.4571420e-03  -1.7666820e-03  -1.1110310e-03  -5.2213700e-04  -2.1458000e-05   3.7860900e-04   6.7424800e-04   8.6879700e-04   9.7131700e-04  -9.9420500e-04  -9.5415100e-04  -8.6641300e-04  -7.4720400e-04  -6.1178200e-04  -4.7254600e-04  -3.3903100e-04  -2.1886800e-04  -1.1634800e-04  -3.4332000e-05   2.7180000e-05   6.9618000e-05   9.5367000e-05   1.0681200e-04   1.0824200e-04   1.0156600e-04  -9.0599000e-05  -7.6771000e-05  -6.2943000e-05  -4.9591000e-05  -3.7670000e-05  -2.7657000e-05  -1.9550000e-05  -1.3828000e-05  -9.0600000e-06  -6.1990000e-06  -3.8150000e-06  -2.3840000e-06  -1.4310000e-06  -9.5400000e-07  -4.7700000e-07   0.0000000e+00];

  %C = Table_analysis_window();
  %C = C(1:4:end);

  %load h127
  %C = h127.*(1-2*(bitand([0:4:511],64)>0));

  xlen = length(X);
  wlen = length(C); % has to be 128

  subbands = 8;  % subband decimation
  hop = subbands;
  tsteps = floor((xlen - wlen + 1)/hop);

  U = zeros(subbands, tsteps);

  VECTORIZED = 1;

  if VECTORIZED

    i = 2*[0:7]+1;
    k = 4-[0:15];
    %k = [0:63];
    %M = cos((repmat(i',1,16).*repmat(k,8,1))*pi/16);
    % Complex filterbank
    M = exp(sqrt(-1)*(repmat(i',1,16).*repmat(k,8,1))*pi/16);

    for t = 1:floor((xlen-wlen+1)/hop)

            % Extract samples and apply window
            Z = X((t-1)*hop+[1:wlen]) .* C;
            % Sum up with stride 8 (16 blocks)
            Y = sum(reshape(Z,16,8),2);

            % Apply modulation to get 8 subbands
            U(:,t) = M*Y;

    end

  else

    % Old (equivalent, nonvectorized) version

    % Calculate the analysis filter bank coefficients
    for i = 0 : 7,
          for k = 0 : 15,
                  %M(i + 1, k + 1) = cos((2 * i + 1) * (k - 4) * pi / 16);
                  % Complex filterbank
                  M(i + 1, k + 1) = exp(sqrt(-1)*((2 * i + 1) * (k - 4) * pi / 16));
          end
    end


    for t = 1:floor((xlen-wlen+1)/hop)

          Z = X((t-1)*hop+[1:wlen]) .* C;

          % Partial calculation: 16 Yi coefficients
          Y = zeros(1, 16);

          for i = 1 : 16,
                  for j = 0 : 7,
                  Y(i) = Y(i) + Z(i + 16 * j);
                  end
          end

          % Calculate the 8 subband samples Si
          for i = 1 : 8
                  for k = 1 : 16
                          U(i,t) = U(i,t) + M(i, k) * Y(k);
                  end
          end

    end

  end

end
