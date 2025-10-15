function pd_power = pd_output(varargin)
p = inputParser();
p.addRequired('laser'); % laser object
p.addRequired('ring_resonator'); % ring resonator object
p.addRequired('beta'); % modulation depth
p.addRequired('t') % time grid
p.addOptional('bessel_terms',10);
p.parse(varargin{:});

l = p.Results.laser;
rr = p.Results.ring_resonator;
beta = p.Results.beta;
bessel_terms = p.Results.bessel_terms;
t = p.Results.t;


n   = -bessel_terms:bessel_terms;
Jn  = besselj(n, beta);
wn  = l.w_c + n.*l.w_m;
Hn  = rr.Through_Ratio(wn);
coeff  = (Jn .* Hn).';
expMat = exp(1j * (wn.' * t));
E_t    = sqrt(l.P_c) * (coeff.' * expMat);
pd_power = abs(E_t).^2;
end