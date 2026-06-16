function [tfmap, freqs] = timefrequencymap(x, times, pnts, trials, min_freq, max_freq, num_freqs)

% define wavelet parameters
freqs = logspace(log10(min_freq),log10(max_freq),num_freqs);
s = logspace(log10(3),log10(10),num_freqs)./(2*pi*freqs);

% definte convolution parameters
n_wavelet            = length(times);
n_data               = pnts*trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = round((n_wavelet-1)/2);

% get FFT of data
% eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
eegfft = fft(reshape(x,1,pnts*trials),n_conv_pow2);

% initialize
tfmap = zeros(num_freqs, pnts); % frequencies X time X trials

baseidx = dsearchn(times',[-1.0 -0.0]');

% loop through frequencies and compute synchronization
for fi=1:num_freqs
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*freqs(fi).*times) .* exp(-times.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv, pnts, trials)).^2,2);
    tfmap(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end