classdef RingResonator
    properties
        SPEED_OF_LIGHT = 3e8
        n_g = 3.76                % group index
        L = 12.1e-3               % ring length (m)
        k = 0.35                  % field coupling coefficient
        t                          % through field coefficient
        Np_to_dB = 20/log(10);    % field dB to Np
        loss_dB_p_cm = 0.4        % POWER loss in dB/cm
        a                          % round-trip FIELD attenuation factor
        lambda = 1550e-9
        FSR_GHz                    % free spectral range
        beta                       % propagation constant vs omega
        FWHM                       % Full width half max
        Finesse                    % Finesse
        Q_factor                   % Quality Factor
        Freq_Linewidth             % Frequency Linewidth
    end 
    methods
        function obj = RingResonator(varargin)
            obj.t = sqrt(1 - obj.k^2);
            
            alpha_p_dB_per_m = obj.loss_dB_p_cm * 100;             % dB/m (power)
            alpha_field_Np_per_m = alpha_p_dB_per_m / obj.Np_to_dB;% Np/m (field)
            obj.a = exp(-alpha_field_Np_per_m * obj.L);
            
            obj.FSR_GHz = (obj.SPEED_OF_LIGHT/(obj.n_g*obj.L))/1e9;
            obj.beta = @(omega) obj.n_g .* omega ./ obj.SPEED_OF_LIGHT; 
            obj.FWHM = obj.lambda^2*(1 - obj.a*obj.t)/(obj.n_g*obj.L*pi*sqrt(obj.a*obj.t));
            obj.Freq_Linewidth = obj.FWHM*obj.SPEED_OF_LIGHT/obj.lambda^2;
            obj.Finesse = pi*sqrt(obj.a*obj.t)/(1 - obj.a*obj.t);
            
        
        end
        
        function E_out_p_E_in = Through_Ratio(obj, omega)
            
            phase = exp(-1j .* obj.beta(omega) .* obj.L);
            n = obj.t - obj.a .* phase;
            d = 1 - (obj.t .* obj.a) .* phase;
            E_out_p_E_in = n ./ d;
        end
        
        function Tfield = TransmissionField(obj, omega)
            % |E_out/E_in|
            Tfield = abs(obj.Through_Ratio(omega));
        end
        
        function Tpower = Transmission(obj, omega)
            % Power transmission |E_out/E_in|^2
            Tpower = abs(obj.Through_Ratio(omega)).^2;
        end
        
        function PlotTransmission(obj,varargin)
            p = inputParser; p.addParameter('fsr_span',10); p.parse(varargin{:});
            fsr_span = p.Results.fsr_span;
            x = linspace(-fsr_span*obj.FSR_GHz, fsr_span*obj.FSR_GHz, 100000)*1e9*2*pi;
            sweep = x;
            Tpower = obj.Transmission(sweep);
            
            figure;
            plot(sweep, Tpower, 'LineWidth', 1.5); grid on;
            xlabel('\Delta\nu (GHz)');
            ylabel('Power transmission  |E_{out}/E_{in}|^2');
            title(sprintf('Ring Through-Port Power vs Detuning (FSR = %.2f GHz)', obj.FSR_GHz));
            legend('Through-port power', 'Location', 'best')
        end
        
        function PlotFieldPhase(obj,varargin)
            p = inputParser; p.addParameter('fsr_span',10); p.parse(varargin{:});
            fsr_span = p.Results.fsr_span;
            x = linspace(-fsr_span*obj.FSR_GHz, fsr_span*obj.FSR_GHz, 100000)*1e9*2*pi;
           
            sweep = x;
            Tfield = obj.Through_Ratio(sweep);
            
            figure;
            plot(sweep, angle(Tfield)*180/pi, 'LineWidth', 1.5); grid on;
            xlabel('\Delta\nu (GHz)');
            ylabel('Phase Degrees');
            title(sprintf('Ring Through-Port Power vs Detuning (FSR = %.2f GHz)', obj.FSR_GHz));
            legend('Through-port power', 'Location', 'best')
        end
        function PlotFieldAndPhase(obj, varargin)
            
            p = inputParser; p.addParameter('fsr_span',10); p.parse(varargin{:});
            fsr_span = p.Results.fsr_span;
            
            x = linspace(-fsr_span*obj.FSR_GHz, fsr_span*obj.FSR_GHz, 100000)*1e9*2*pi;
            %omega_0 =  2*pi*obj.SPEED_OF_LIGHT/obj.lambda;
            sweep = x;
            
            
            E = obj.Through_Ratio(sweep);
            mag = abs(E);
            ph_deg = angle(E)*180/pi;                  
            ph_deg = mod(ph_deg + 180, 360) - 180;      % wrap to (-180, 180]
     
            ph_deg(mag < 1e-8) = NaN;
            
            figure; grid on; hold on;
            yyaxis left
            plot(sweep, mag.^2, 'LineWidth', 1.5);
            ylabel('|E_{out}/E_{in}|');
            
            yyaxis right
            plot(sweep, ph_deg, '--', 'LineWidth', 1.25);
            ylabel('Phase (deg)');
            
            xlabel('\Delta\nu (GHz)');
            title(sprintf('Through-Port Field & Phase vs Detuning (FSR = %.2f GHz)', obj.FSR_GHz));
            legend({'Field magnitude','Phase'}, 'Location','best');
        end
    end
end
