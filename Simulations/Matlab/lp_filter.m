function filtered_signal = lp_filter(varargin)
p = inputParser;
p.addRequired('signal');
p.addRequired('laser') % laser
p.addRequired('t'); % time grid
p.parse(varargin{:});

signal = p.Results.signal; 
l = p.Results.laser;
t = p.Results.t;

fc_lpf = 0.2*l.f_m; % make the cutoff 20% on the modulation frequency

fs = 1/mean(diff(t)); % sample rate
Wn = fc_lpf/(fs/2);
[b,a]  = butter(2, Wn); % 2nd order butterworth
filtered_signal   = filtfilt(b,a,signal);
end