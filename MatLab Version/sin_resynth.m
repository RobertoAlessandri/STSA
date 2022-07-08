function s=sin_resynth(sinan_freqs, sinan_amps, sinan_phs, Fs, SpF)
%%% define controls (analysis matrices sinan_freqs and
%%% sinan_amps have been created in the analysis phase)
npart=size(sinan_amps,1); %number of analyzed partials 
t0=0; %initial time
%%% compute sound %%%
s=0; %initialize output signal
for i=1:npart %generate all partials and sum
    s=s+sinosc(t0,sinan_amps(i,:),sinan_freqs(i,:),sinan_phs(i,:), Fs, SpF);
end
