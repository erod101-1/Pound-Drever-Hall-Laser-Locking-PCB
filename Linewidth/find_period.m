function f_Hz = find_period(S,psd_dbm_p_hz)
%FIND_PERIOD Estimate the dominant period (in Hz) of a PSD-vs-frequency trace via autocorrelation.
%
%   f_Hz = FIND_PERIOD(S, psd_dbm_p_hz)
%
%   Inputs
%   ------
%   S : struct
%       Must contain field:
%           S.TraceDataFrequency  (Nx1 or 1xN double) frequency axis in Hz, evenly (or near-evenly) spaced.
%
%   psd_dbm_p_hz : (Nx1 or 1xN double)
%       PSD values (e.g., dBm/Hz) sampled on the same grid as S.TraceDataFrequency.
%
%   Output
%   ------
%   f_Hz : double
%       Estimated period expressed in Hz (e.g., FSR-like spacing), computed as:
%           period_samples * mean(diff(S.TraceDataFrequency)).
%
%   Notes
%   -----
%   - The signal is detrended by subtracting a moving average (movmean) before autocorrelation.
%   - Autocorrelation is computed with xcorr(...,'coeff'), then only positive lags are searched.
%   - The first detected peak in the autocorrelation is used as the period in samples.
%
%   Example
%   -------
%   f_Hz = find_period(S, psd_dbm_p_hz);
%% Remove the sincs 
x = psd_dbm_p_hz(:);

% subtract av
W = 200;  % <-- tune this (often ~1 period in samples, or larger to remove slow drift)
x = x - movmean(x, W, "Endpoints","shrink");
x = x - mean(x);

% compute autocorelation
[ac,lags] = xcorr(x, x, "coeff");

% keep the positive peaks
pos = lags > 0;
acp = ac(pos);
lagp = lags(pos);

%% Plot Autocorrelation
figure;
plot(lagp, acp);
xlabel('Units (arb)');
ylabel('Units (arb)');
title('Autocorelation');

% find the peak near the period
[~,locs] = findpeaks(acp, lagp);
f_samp = locs(1);                
df = mean(diff(S.TraceDataFrequency)); % frequency step from data;
f_Hz = f_samp*df;
end