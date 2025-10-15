classdef Laser
    properties
        SPEED_OF_LIGHT = 3e8
        lambda double % carrier wavelength
        f_c double % carrier frequency
        w_c double
        f_m double % modulation frquency
        w_m double
        P_c double % carrier power
        % signal
    end
    methods

        function obj = Laser(varargin)
            p = inputParser;
            p.addParameter('lambda',1550e-9);
            p.addParameter('P_c',5e-3);
            p.addParameter('f_m',1e9);
            p.parse(varargin{:});
            obj.lambda = p.Results.lambda;
            obj.P_c = p.Results.P_c; % carrier power W
            obj.f_m = p.Results.f_m; % modulation frequency Hz
            obj.f_c = obj.SPEED_OF_LIGHT/obj.lambda; % carrier frequency Hz
            obj.w_m = 2*pi*obj.f_m;
            obj.w_c = 2*pi*obj.f_c;

        end
        function [spectrumArr,spectrumAmp] = sidebands(obj,terms,beta)
            spectrumArr = -terms:terms;
            ampArr = sqrt(obj.P_c)*besselj(spectrumArr, beta);
            spectrumAmp = abs(ampArr).^2;
            spectrumArr = obj.w_c - spectrumArr*obj.w_m;

            % Plot spectrum
            f_k = spectrumArr/(2*pi);
            figure; stem(f_k/1e9, spectrumAmp*1000, 'filled'); grid on
            xlabel('Frequency (GHz)'); ylabel('Line power (mW)');
            title('PDH line spectrum');
        end
    end
end