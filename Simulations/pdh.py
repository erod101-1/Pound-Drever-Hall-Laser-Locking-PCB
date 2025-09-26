from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 3e8 # speed of light

@dataclass
class laser:
    P_c : float
    f_c_THz : float # THz
    f_c_Hz : float = field(init=False)
    l_0: float  = field(init=False) 
    omega_c: float = field(init=False)  # rad/s
    def __post_init__(self):
        self.f_c_Hz = self.f_c_THz * 1e12 # convert THz to Hz
        self.l_0 = SPEED_OF_LIGHT / self.f_c_Hz 
        self.omega_c = 2*np.pi * self.f_c_Hz # 2*pi*Hz
    
@dataclass
class eom:
    f_m : float # modulation frequency
    omega_m : float = field(init=False) # sideband angular frequency
    P_s : float # side band power
    def __post_init__(self):
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
    
    # get reflection coeffient (from black paper)
    def reflection_coefficient(self,omega : np.linspace) -> float:
        n = self.r*(np.exp(1j*omega/self.fsr)-1)
        d = 1 - self.r*self.r*np.exp(1j*omega/self.fsr)
        return n/d

# Laser... -
mcavity = cavity(r=0.95, L=0.02)
mlaser = laser(P_c=1e-3,f_c_THz=193.15)
meom = eom(f_m=0.1e9, P_s=0.5e-3)

omega_space = np.linspace(mlaser.omega_c-3*meom.omega_m, \
                                       mlaser.omega_c+3*meom.omega_m, \
                                       int(1e6))


def reflected_power(t):
    # Detuning grid relative to the carrier
    delta = omega_space - mlaser.omega_c
    # Reflection coefficients at carrier and sidebands
    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + meom.omega_m)  # Omega + Omega
    Fm = mcavity.reflection_coefficient(delta - meom.omega_m)  # Omega − Omega
    Pc = mlaser.P_c
    Ps = meom.P_s
    Omega  = meom.omega_m
    # refer to black paper 
    t1 = Pc * (np.abs(F0) ** 2)
    t2 = Ps * ((np.abs(Fp) ** 2) + (np.abs(Fm) ** 2))
    A = F0 * np.conj(Fp) - np.conj(F0) * Fm
    t3 = 2 * np.sqrt(Pc * Ps) * (np.real(A) * np.cos(Omega * t) + np.imag(A) * np.sin(Omega * t))
    return t1 + t2 + t3

def time_averaged_power() -> float:
    delta = omega_space - mlaser.omega_c
    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + meom.omega_m)
    Fm = mcavity.reflection_coefficient(delta - meom.omega_m)
    return (mlaser.P_c*np.abs(F0)**2 + meom.P_s*(np.abs(Fp)**2 + np.abs(Fm)**2))


def error_signal(demod_phase=np.pi/2):

    delta = omega_space - mlaser.omega_c
    Om = meom.omega_m

    F0 = mcavity.reflection_coefficient(delta)
    Fp = mcavity.reflection_coefficient(delta + Om)
    Fm = mcavity.reflection_coefficient(delta - Om)

    Pc, Ps = mlaser.P_c, meom.P_s

    A = F0*np.conj(Fp) - np.conj(F0)*Fm

    Verr = 2*np.sqrt(Pc*Ps) * (
        np.real(A)*np.cos(demod_phase) + np.imag(A)*np.sin(demod_phase)
    )
    return delta, Verr


delta, Verr = error_signal()

plt.figure(figsize=(9,6))
plt.subplot(2,1,1)
plt.plot(delta, time_averaged_power())
plt.title("Reflected Power (time-averaged) vs Detuning")
plt.ylabel("Power (W)")
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(delta, Verr)
plt.title("PDH Error Signal vs Detuning")
plt.xlabel("Detuning (rad/s)")
plt.ylabel("Error (arb. units)")
plt.grid(True)
plt.tight_layout()
plt.show()


