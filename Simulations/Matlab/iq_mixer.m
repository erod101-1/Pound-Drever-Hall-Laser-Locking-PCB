function [I_raw,Q_raw] = iq_mixer(varargin)
    p = inputParser;
    p.addRequired('pd_output'); % photodiode output signal
    p.addRequired('laser'); % object contains modualtion frequency
    p.addRequired('t'); % time grid
    p.addOptional('phi',0);
    p.parse(varargin{:});
    pd_output = p.Results.pd_output;
    l = p.Results.laser;
    t = p.Results.t;
    phi = p.Results.phi;
    LOc = cos(l.w_m*t - phi); % I
    LOs = sin(l.w_m*t - phi); % Q
    I_raw = real(pd_output .* LOc);
    Q_raw = real(pd_output .* LOs);
end