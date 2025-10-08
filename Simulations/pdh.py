from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 3e8 # speed of light
DBM_TO_MW = lambda x : 10
@dataclass
class laser:
    P_c : float # Carrier Power
    f_c_THz : float # THz
    f_m : float # modulation frequency
    omega_m : float = field(init=False) # angular modualtion frequency
    f_c_Hz : float = field(init=False) # carrier power in Hz
    P_s : float # sideband power in MW
    l_0: float  = field(init=False)  # wavelength of the carrier
    omega_c: float = field(init=False)  # angular carrier frequency
    def __post_init__(self):
        self.f_c_Hz = self.f_c_THz * 1e12 # convert THz to Hz
        self.l_0 = SPEED_OF_LIGHT / self.f_c_Hz 
        self.omega_c = 2*np.pi * self.f_c_Hz # 2*pi*Hz
        self.omega_m = 2*np.pi*self.f_m 

    
@dataclass 
class cavity:
    r : float  # reflection coefficient
    L : float  # length of the cavity
    fsr : float = field(init=False) # free spectral range
    def __post_init__(self):
        self.fsr = SPEED_OF_LIGHT/(2*self.L)

    def reflection_coefficient(self,omega : np.linspace) -> float:
        n = self.r*(np.exp(1j*omega/self.fsr)-1)
        d = 1 - self.r*self.r*np.exp(1j*omega/self.fsr)
        return n/d
    
# Laser... -
mcavity = cavity(r=0.95, L=0.02)
mlaser = laser(P_c=1e-3,f_c_THz=193.15,f_m=0.1e9, P_s=0.5e-3)
omega_space = np.linspace(mlaser.omega_c-3*mlaser.omega_m, \
                                       mlaser.omega_c+3*mlaser.omega_m, \
                                       int(1e6))



def error_signal(_laser: laser,
                 omega_grid: np.ndarray,
                 demod_phase: float,
                 phase_offset_rad: float = 0.0,
                 path_delay_s: float = 0.0):
    
    delta = omega_grid - _laser.omega_c
    Om    = _laser.omega_m

    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + Om)
    Fm = mcavity.reflection_coefficient(delta - Om)

    Pc, Ps = _laser.P_c, _laser.P_s

    A = F0 * np.conj(Fp) - np.conj(F0) * Fm

    # effective phase = demod phase + internal fixed offset + delay-induced phase
    phi_eff = demod_phase + phase_offset_rad + (Om * path_delay_s)

    Verr = 2.0 * np.sqrt(Pc * Ps) * (
        np.real(A) * np.cos(phi_eff) + np.imag(A) * np.sin(phi_eff)
    )

    # ensure strictly real (numeric safety)
    return delta, np.asarray(np.real_if_close(Verr), dtype=float)

def get_error_slope_magnitude(error_signal : np.ndarray):
    error_signal_window = error_signal[int(len(error_signal)/2) : int(len(error_signal)/2)+10]
    dx = 1/int(1e6)
    return np.average(np.gradient(error_signal_window,dx))

def plot_error_signal(Verr,delta):
    plt.plot(delta, Verr)
    plt.title("PDH Error Signal vs Detuning")
    plt.xlabel("Detuning (rad/s)")
    plt.ylabel("Error (arb. units)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

delta, Verr = error_signal(_laser=mlaser,omega_grid=omega_space,demod_phase=0,phase_offset_rad=0*np.pi/180)
plot_error_signal(Verr,delta)

#max_error_magnitude = 0
#for phase in np.linspace(0,180,180):
#    delta, Verr = error_signal(_laser=mlaser,omega_grid=omega_space,demod_phase=0,phase_offset_rad=phase*np.pi/180)
#    
#    if get_error_slope_magnitude(error_signal=Verr) > max_error_magnitude:
#        max_error_magnitude = get_error_slope_magnitude(error_signal=Verr)
#        max_phase = phase
#
#print(max_error_magnitude)
#print(max_phase)








'''
Dont really care about power
def reflected_power(t):
    # Detuning grid relative to the carrier
    delta = omega_space - mlaser.omega_c
    # Reflection coefficients at carrier and sidebands
    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + mlaser.omega_m)  # Omega + Omega
    Fm = mcavity.reflection_coefficient(delta - mlaser.omega_m)  # Omega − Omega
    Pc = mlaser.P_c
    Ps = mlaser.P_s
    Omega  = mlaser.omega_m
    # refer to black paper 
    t1 = Pc * (np.abs(F0) ** 2)
    t2 = Ps * ((np.abs(Fp) ** 2) + (np.abs(Fm) ** 2))
    A = F0 * np.conj(Fp) - np.conj(F0) * Fm
    t3 = 2 * np.sqrt(Pc * Ps) * (np.real(A) * np.cos(Omega * t) + np.imag(A) * np.sin(Omega * t))
    return t1 + t2 + t3

def time_averaged_power() -> float:
    delta = omega_space - mlaser.omega_c
    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + mlaser.omega_m)
    Fm = mcavity.reflection_coefficient(delta - mlaser.omega_m)
    return (mlaser.P_c*np.abs(F0)**2 + mlaser.P_s*(np.abs(Fp)**2 + np.abs(Fm)**2))
'''