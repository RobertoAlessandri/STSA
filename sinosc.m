function s = sinosc(t0, a, f, ph0, Fs, SpF)

% φ(t) = φ(tk) + 2πf0(tk)(t − tk) + 2πFc[f0(tk+1) − f0(tk)] * ((t - tk)**2) / 2 

% Fs sample rate
% SpF samples-per-frame 
nframes = length(a); % total number of frames
if (length(f)==1)
    f = f*ones(1,nframes);
end
if (length(f) ~= nframes)
    error('wrong f length!')
end
s = zeros(1, nframes*SpF); % initialize signal vector to 0
lasta = a(1); lastf = f(1); lastph = ph0; % initialize amplitude,frequency, phase

for i = 1:nframes % cycle on the frames
    naux = 1:SpF; % count samples within frame
    % Compute amplitudes and phases within frame %
    ampl = lasta + (a(i) - lasta)/SpF.*naux;
    phase = lastph + pi/Fs.*naux.*(2*lastf + (1/SpF)*(f(i)-lastf).*naux);
    % read from table %
    s(((i-1)*SpF+1):i*SpF) = ampl.*cos(phase); % read from table
    % save last values of amplitude, frequenc, phase %
    lasta = a(i); lastf = f(i);
    lastph = phase(SpF);
end
s = [zeros(1, round(t0*Fs)) s]; % add initial silence of t0 sec.
end

