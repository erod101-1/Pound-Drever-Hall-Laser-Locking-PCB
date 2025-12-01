function [f_Hz, S_nu_f] = get_noise_psd(S)
%GET_NOISE_PSD Estimate frequency-noise PSD from SA trace, remove sinc imprint, and plot beta line.
% Outputs:
%   f_Hz   : frequency axis (Hz)
%   S_nu_f : frequency-noise PSD (Hz^2/Hz), linear

%% Frequency axis (assume S.TraceDataFrequency is in Hz)
f_Hz = S.TraceDataFrequency(:);
f_MHz = f_Hz/1e6;

%% Plot the original data
figure;
plot(f_MHz, S.TraceDataPower);
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dBm/RBW)');
title('Original Spectrum Analyzer Data');

%% Normalize data by the RBW (to dBm/Hz)
S.TraceDataPower = S.TraceDataPower(:);
S.TraceDataPower = [S.TraceDataPower(2); S.TraceDataPower(2:end)];

RBW_Hz = S.SetupInfo.TRACEA.RBW * 10e3;                % (kept as your original)
psd_dbm_p_hz = S.TraceDataPower - 10*log10(RBW_Hz);

figure;
plot(f_MHz, psd_dbm_p_hz);
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dBm/Hz)');
title('Spectrum Analyzer Data RBW Normalized');

%% Remove the sincs...
f_period_Hz = find_period(S, psd_dbm_p_hz);
T_period_s = 1/f_period_Hz;

I0 = (5e-3)/2;
K0 = 1
K1 = (2*pi*T_period_s)^2;

Sy = 1e-3 * 10.^(psd_dbm_p_hz/10);                     % W/Hz (linear)
H  = sinc(f_Hz * T_period_s);                          % MATLAB sinc = sin(pi x)/(pi x)
H2 = H.^2;

floorMag = 2e-2;                                       % tune
den = max(H2, floorMag^2);

S_nu_f = Sy ./ (K0 .* K1 .* den);                      % Hz^2/Hz (linear)
S_nu_dB = 10*log10(S_nu_f);

%% Plot S_nu(f) with beta line
mask = isfinite(S_nu_f) & (f_Hz > 0) & (S_nu_f > 0);

beta_line = (8*log(2)/pi^2) * f_Hz;                    % Hz^2/Hz

figure; hold on; grid on;
loglog(f_Hz(mask), S_nu_f(mask));
loglog(f_Hz(mask), beta_line(mask), '--');
xlabel('Frequency offset f (Hz)');
ylabel('S_\nu(f) (Hz^2/Hz)');
title('Frequency Noise PSD with \beta-line');
legend('S_\nu(f)', '\beta-line: (8 ln2/\pi^2) f', 'Location','best');
hold off;

%% Optional: same plot in dB (if you like)
figure; hold on; grid on;
semilogx(f_Hz(mask), 10*log10(S_nu_f(mask)));
semilogx(f_Hz(mask), 10*log10(beta_line(mask)), '--');
xlabel('Frequency offset f (Hz)');
ylabel('Level (dB[Hz^2/Hz])');
title('Frequency Noise PSD (dB) with \beta-line');
legend('S_\nu(f)', '\beta-line', 'Location','best');
hold off;

end
