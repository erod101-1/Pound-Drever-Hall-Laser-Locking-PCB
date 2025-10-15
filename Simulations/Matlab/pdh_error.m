%% Inputs
MODULATION_FREQUENCY = 1e9;
bessel_terms = 10;
beta = 1; % modulation depth
td_periods = 10; % time domain periods to show
tm = 1/MODULATION_FREQUENCY;
samples_per_period = 4096;
phi = 0; % delay
fc_lpf = 0.2*l.f_m;
%% Setup
% cavity
rr = RingResonator();

% find the nearest resonance frequency
target_lambda_m = 1550e-9;
m = (rr.SPEED_OF_LIGHT/target_lambda_m)/(rr.FSR_GHz*1e9);
m = round(m); % round to the nearest decimal
closest_resonant_frequency = m*rr.FSR_GHz*1e9;
fprintf('Closest Resonant Frequency: %f THz\n',closest_resonant_frequency/1e12);
omega_closest = closest_resonant_frequency*2*pi;
lambda_closest = rr.SPEED_OF_LIGHT/closest_resonant_frequency;

% instaniate laser
l = Laser('lambda',lambda_closest,'f_m',MODULATION_FREQUENCY);
% time domain grid
t = linspace(0,td_periods*tm,td_periods*samples_per_period);
fs = 1/mean(diff(t)); % sample rate

%% Photodiode response
pd_power = pd_output(l,rr,beta,t);


%% IQ Mixer
[I_raw,Q_raw] = iq_mixer(pd_power,l,t,'phi',0);


IQ_Raw = I_raw + Q_raw;

%% Filter
I_filtered = lp_filter(I_raw,l,t);
Q_filtered = lp_filter(Q_raw,l,t);
opt_angle = atan2(mean(Q_filtered),mean(I_filtered))*180/pi;

fprintf('Optimum demod angle: %f',opt_angle);

% pd power
figure('Name','Time-Domain PD Response');
subplot(4,1,1);
plot(t*1e9, pd_power, 'LineWidth', 1.3); grid on;
xlabel('Time (ns)'); ylabel('PD power  \propto |E(t)|^2');
title(sprintf('PD Power (\\beta=%.2f, terms=%d, f_m=%.3f GHz)', ...
    beta, bessel_terms, l.f_m/1e9));


% IQ
subplot(4,1,2);
xlabel('Time (ns)'); ylabel('I and Q raw');
plot(t*1e9, I_raw, 'LineWidth', 1.2); hold on;
plot(t*1e9, Q_raw, 'LineWidth', 1.2); grid on;
legend('I_{raw}','Q_{raw}','Location','best');
title('Raw IQ Mixer Outputs at f_m');

% filter
subplot(4,1,3);
plot(t*1e9, I_filtered, 'LineWidth', 1.4); hold on;
plot(t*1e9, Q_filtered, 'LineWidth', 1.4); hold on;
grid on; xlabel('Time (ns)'); ylabel('Baseband');
legend('(I)_{LPF}','(Q)_{LPF}','Location','best');
title(sprintf('LPF Baseband (Butterworth f_c=%.2f MHz)', fc_lpf/1e6));

% opt_angle*pi/180
demod_opt =  opt_angle*pi/180;
IQ_bb = I_filtered*cos(demod_opt) + Q_filtered*sin(demod_opt);
subplot(4,1,4);
plot(t*1e9, IQ_bb, 'LineWidth', 1.4); hold on;
grid on; xlabel('Time (ns)'); ylabel('Baseband');
legend('(IQ)_{BB}','Location','best');
title(sprintf('I + Q at optimal demod f_c=%.2f MHz', fc_lpf/1e6));



%
%%% -------- Frequency Sweeps — sweep beta, recompute optimal demod phase per beta ----
%Ndnu   = 500;
%dnu_Hz = linspace(-rr.Freq_Linewidth, rr.Freq_Linewidth, Ndnu);  % Hz detuning (already in Hz)
%domega = 2*pi*dnu_Hz;                                            % rad/s detuning
%
%beta_list = [0.1 0.3 0.7 1.0 2.0 3.0];
%
%error_by_beta = zeros(numel(beta_list), Ndnu);
%phi_opt_deg   = zeros(1, numel(beta_list));
%
%% index of detuning closest to zero (in rad/s)
%[~, k0] = min(abs(domega));
%Kfit = 3;                                % half-window for slope fit (use k0-Kfit : k0+Kfit)
%
%for ib = 1:numel(beta_list)
%    beta_k = beta_list(ib);
%    fprintf('\n BETA STEP: %d \n\t',ib);
%    % --- DC I/Q vs detuning for this beta ---
%    DC_I = zeros(1, Ndnu);
%    DC_Q = zeros(1, Ndnu);
%
%    for k = 1:Ndnu
%        fprintf('DNU STEP: %d \n',k);
%        % Retune laser by domega(k)
%        l_k        = l;
%        l_k.w_c    = l.w_c + domega(k);
%        l_k.f_c    = l_k.w_c/(2*pi);
%        l_k.lambda = l_k.SPEED_OF_LIGHT / l_k.f_c;
%
%        % Photodiode @ this detuning and beta
%        pd_k = pd_output(l_k, rr, beta_k, t);
%
%        % IQ mix (LO phase = 0)
%        [I_raw_k, Q_raw_k] = iq_mixer(pd_k, l_k, t, 'phi', 0);
%
%        % Low-pass, then DC
%        I_bb_k = lp_filter(I_raw_k, l_k, t);
%        Q_bb_k = lp_filter(Q_raw_k, l_k, t);
%        DC_I(k) = mean(I_bb_k);
%        DC_Q(k) = mean(Q_bb_k);
%    end
%
%    % --- Optimal demod phase at center (max-slope), using slopes vs domega ---
%    i1 = max(1, k0-Kfit);  i2 = min(Ndnu, k0+Kfit);
%    x  = domega(i1:i2);                 % rad/s
%    pI = polyfit(x, DC_I(i1:i2), 1);    % [mI_rad, bI]
%    pQ = polyfit(x, DC_Q(i1:i2), 1);    % [mQ_rad, bQ]
%    mI_rad = pI(1);  mQ_rad = pQ(1);
%
%    phi_opt = atan2(mQ_rad, mI_rad);    % radians
%    phi_opt_deg(ib) = rad2deg(phi_opt);
%
%    % --- Discriminator for all detunings at this fixed φ ---
%    error_by_beta(ib, :) = DC_I.*cos(phi_opt) + DC_Q.*sin(phi_opt);
%end
%
%% ---- Plot all beta curves on one figure (x = Δω in rad/s) ----
%figure('Name','PDH Error vs Detuning — Beta Sweep (φ re-estimated per β)'); hold on; grid on;
%cols = lines(numel(beta_list));
%for ib = 1:numel(beta_list)
%    plot(domega, error_by_beta(ib,:), 'LineWidth', 1.8, 'Color', cols(ib,:));
%end
%xlabel('\Delta\omega (rad/s)');
%ylabel('Error Signal (DC)');
%title(sprintf('PDH Error vs Detuning (f_m = %.3f GHz, \\phi_{opt} per \\beta)', l.f_m/1e9));
%leg = arrayfun(@(b,ph) sprintf('\\beta = %.2f  (\\phi_{opt}=%.1f^\\circ)', b, ph), ...
%               beta_list, phi_opt_deg, 'UniformOutput', false);
%legend(leg, 'Location', 'best'); hold off;
