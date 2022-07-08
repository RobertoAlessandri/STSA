clearvars;
close all;
clc;

% Load clean speech signal
[sig_m, Fs] = audioread('hiss_fast_cle.wav');
sig_m = sig_m(:, 1);

% Remove the mean to be sure to have a zero-mean signal
mean_sig = mean(sig_m);
sig = sig_m - mean_sig;

L = length(sig); % Signal length
pow2 = floor(log2(L));
L_nfft = 2^pow2;

%% Short-time processing params
winlen = 2048; % Length of the analysis window
win = hann(winlen, 'periodic'); % Window
hopsize = winlen / 4;
numHopsPerFrame = winlen / hopsize; 
n_frames = floor((L-winlen)/hopsize) + 1; % number of time frames
nfft = winlen; % Number of FFT frequency bins

figure(1)
stft(sig_m, Fs, 'FrequencyRange', 'onesided', 'FFTLength', nfft)

%% Noise Estimation
% We search for the more suited half-second frame to use for the noise
% estimatation.
% We search inside the audio signal for the frame with the best combination
% of zero-crossing rate (high) and short-time energy.
only_noise_samples = Fs/2;

hopsize_noise = only_noise_samples;
n_frames_noise = floor((L-only_noise_samples)/hopsize_noise) + 1; % number of time frames

% These vectors are used for a later representation of the signal
zcr_vector = zeros(1, n_frames_noise);
ste_vector = zeros(1, n_frames_noise);

% We use as starting values for the later searching of the optimum ones, 
% the average zcr and ste of the whole audio signal
zcr_whole = (1/(L-1))*sum(abs(diff(sign(sig))));
ste_whole = sum(sig.^2) / L;
zcr = zcr_whole;
ste = ste_whole;
noise_frame = 0; % we are going to save the frame in order to have a better 
                 % understanding if the algorithm pointed out a reasonable
                 % frame

for nn=0:n_frames_noise-1
    % Segment the corrupted signal x 
    x_n = sig(nn*hopsize_noise + 1 : nn*hopsize_noise + only_noise_samples);
    
    zcr_frame = (1/(only_noise_samples-1))*sum(abs(diff(sign(x_n))));
    ste_frame = sum(x_n.^2) / only_noise_samples;
    zcr_vector(nn + 1) = zcr_frame;
    ste_vector(nn + 1) = ste_frame;
    % We assume that noisiest part of an audio signal are those with low
    % ste and high zcr
    if (zcr_frame > zcr && ste_frame < ste)
        noise_frame = nn;
        zcr = zcr_frame;
        ste = ste_frame;
        x_noise = x_n;
    end
end

% Vector normalization just for visualization
zcr_vector_n = zcr_vector / max(zcr_vector);
ste_vector_n = ste_vector / max(ste_vector);
t_nn = 0:n_frames_noise - 1;

figure(2); 
plot(t_nn, ste_vector_n');
hold on;
plot(t_nn, zcr_vector_n', '--');
title('ZCR and STE raw evolution over time')
xlabel('Time [2 frames per s]');
ylabel('ZCR and STE normalized');
legend('ste','zcr')
disp(['From figure 1 you can see that the noisiest frame is number: ' num2str(noise_frame)]);

% Noise Power Spectrum Estimate 
r_v = xcorr(x_noise);
P_v = abs(fft(r_v, nfft));
P_v_ph = angle(fft(r_v, nfft));
mean_noise = mean(sig_m(noise_frame*hopsize_noise + 1 : noise_frame*hopsize_noise + only_noise_samples));

%% Suppresion Rule Parameters

alpha = [1, 2, 3, 4, 4.5, 6, 8];
beta = [0.1, 0.5, 0.75, 1, 2];
a = 1;
b = 1;
c = 1;

frames_ratio = n_frames / n_frames_noise;

out = zeros(L, 1);
out_resynth = zeros(L, 1);
sinan_freqs = zeros(winlen, 1);
sinan_amps = zeros(winlen, 1);
sinan_phs = zeros(winlen, 1);
G_min = 1e-8;
for i = 1: winlen - 1
    sinan_freqs(i,1) = i * Fs / winlen;
end    
for n=0:n_frames-1
    % Segment the corrupted signal x 
    x_n = sig(n*hopsize + 1 : n*hopsize + winlen);
    x_n_w = win .* x_n; % Windowing
    X_n = fft(x_n_w, nfft); % FFT
    P_X = abs(X_n).^2; % Power Spectrum Estimate
    P_X_min = min(P_X(P_X>0)) / 10;
    P_X(P_X == 0) =  P_X_min; % 1e-50 avoids numerical issues when P_x==0
    N = length(X_n);

    % Gain Estimation
    Z_B = (P_X) - a * (P_v.^c);
    Z_B(Z_B<0) = 0; % This should be considered as the Power Spectrum of the "clean "
                    % signal, so as Power Spectrum, it can't
                    % have negative values (which may appear due to random
                    % fluctuations)
    Z = Z_B ./ (P_X);
    G = (Z).^b;
    %sinan_amps = abs(G);
    %sinan_phs = angle(G);
    if (min(G) > 0)
        G_min = min(G(G>0))/10;
    end    
    disp(['G_min: ' G_min]);
    G(G==0) = G_min;
    G = 1./G;
    G = G/(max(G));
    
    
    % Apply the Gain to the frame spectrum
    S = X_n .* (0.1 * G);
    sinan_amps = abs(S);
    sinan_phs = angle(S);
    % Compute the inverse FFT of the filtered frame. This corresponds to an
    % estimate of the clean speech.
    s_n = ifft(S, nfft);
    s_n_resynth = sin_resynth(sinan_freqs, sinan_amps, sinan_phs, Fs, winlen)';

    % Window the filtered time-domain frame
    s_n_w = win .* s_n;
    s_n_resynth_w = win .* s_n_resynth;
    

    % OLA
    out(n*hopsize + 1 : n*hopsize + winlen) = out(n*hopsize + 1 : n*hopsize + winlen) + s_n_w;
    out_resynth(n*hopsize + 1 : n*hopsize + winlen) = out_resynth(n*hopsize + 1 : n*hopsize + winlen) +  s_n_resynth_w;
    % Save the filter to visualize how it attenuates what should be the
    % noisiest frame
    if (n == floor(noise_frame*frames_ratio))
        clean_frame = n;
        G_mid = G;
        [G__, omega] = freqz(Z_B, P_X, N);%, 'whole');
    end

end

%% Final Results


MSE = sum((out - sig).^2)/ length((out - sig).^2);
var_out = var(out);
var_in = var(sig);

% Power Spectrum of original and cleaned audio signal
P_S_total = abs(fft(out, L_nfft)).^2;
P_X_total = abs(fft(sig, L_nfft)).^2;

% Power Spectrum Noise Estimet for the clean signal
v_unclean = out(noise_frame*hopsize_noise + 1 : noise_frame*hopsize_noise + only_noise_samples);
r_v_unclean = xcorr(v_unclean);
P_v_new = abs(fft(r_v, L_nfft)) + 1e-50;
P_v_unclean_new = abs(fft(r_v_unclean, L_nfft)) + 1e-50;

% Signal-to-Noise Ratio before and after attenuation
SNR_in = pow2db(mean((P_X_total ./ P_v_new)));
SNR_out = pow2db(mean((P_S_total ./ P_v_unclean_new)));
disp(['SNR_in: ' num2str(SNR_in)]);
disp(['SNR_out: ' num2str(SNR_out)]);


out = out + mean_sig + mean_noise;
out_resynth = out_resynth + mean_sig + mean_noise;

real_out = out + out_resynth;

figure(3)
stft(out, Fs, 'FrequencyRange', 'onesided', 'FFTLength', nfft)
%[sinan_amps, sinan_freqs] = stft(out, Fs, 'FrequencyRange', 'onesided', 'FFTLength', nfft);
%out_resynth = sin_resynth(sinan_freqs, sinan_amps);


soundsc(real_out, Fs);

omega = (omega / pi) * Fs;

figure(4); 
plot(omega, G_mid);
title('Filter at the noisiest frame')
xlabel('Frequency [Hz]');
ylabel('Power Spectrum normalized');
legend('Suppresion Rule')