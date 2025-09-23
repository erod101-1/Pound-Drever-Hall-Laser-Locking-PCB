from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt


class CONSTANTS:
    SPEED_OF_LIGHT = 3e8 # speed of light

@dataclass
class laser(CONSTANTS):
    l_0: float
    P_c: float
    f_c: float = field(init=False)      # Hz
    omega_c: float = field(init=False)  # rad/s
    def __post_init__(self):
        self.f_c = self.SPEED_OF_LIGHT / self.l_0
        self.omega_c = 2*np.pi * self.f_c
        

@dataclass
class eom(CONSTANTS):
    f_m : float # modulation frequency
    omega_m : float = field(init=False) # sideband angular frequency
    P_s : float # side band power
    def __post_init__(self):
        self.omega_m = 2*np.pi*self.f_m
@dataclass 
class cavity(CONSTANTS):
    r : float  # reflection coefficient
    L : float  # length of the cavity
    fsr : float = field(init=False) # free spectral range
    def __post_init__(self):
        self.fsr = self.SPEED_OF_LIGHT/(2*self.L)

class ReflectionCoeff(cavity):
    def __init__(self,**kwargs):
        c =kwargs.get('cavity',cavity(r=0.95, L=0.02))
        super().__init__(r=c.r,L=c.L)

    # get reflection coeffient (from black paper)
    def reflection_coefficient(self,omega : np.linspace) -> float:
        n = self.r*(np.exp(1j*omega/self.fsr)-1)
        d = 1 - self.r*self.r*np.exp(1j*omega/self.fsr)
        return n/d

class ReflectedPower(ReflectionCoeff):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        
        self.laser = kwargs.get('laser',laser(l_0=1550e-9, P_c=5e-3))
        self.eom = kwargs.get('eom',eom(f_m=20e6, P_s=0.5e-3))
        self.omega_space = np.linspace(self.laser.omega_c-3*self.eom.omega_m, \
                                       self.laser.omega_c+3*self.eom.omega_m, \
                                       int(1e6))
        
    def reflected_power(self):
        # Detuning grid relative to the carrier
        delta = self.omega_space - self.laser.omega_c

        # Reflection coefficients at carrier and sidebands
        F0 = self.reflection_coefficient(delta)
        Fp = self.reflection_coefficient(delta + self.eom.omega_m)  # Omega + Omega
        Fm = self.reflection_coefficient(delta - self.eom.omega_m)  # Omega − Omega

        Pc = self.laser.P_c
        Ps = self.eom.P_s
        Omega  = self.eom.omega_m

        # refer to black paper 
        t1 = Pc * (np.abs(F0) ** 2)
        t2 = Ps * ((np.abs(Fp) ** 2) + (np.abs(Fm) ** 2))

        A = F0 * np.conj(Fp) - np.conj(F0) * Fm
        t = 0.0
        t3 = 2 * np.sqrt(Pc * Ps) * (np.real(A) * np.cos(Omega * t) + np.imag(A) * np.sin(Omega * t))

        return t1 + t2 + t3

    def error_signal(self, demod_phase=np.pi/2):
        # symmetric detuning grid that includes 0 exactly
        delta = self.omega_space - self.laser.omega_c

        # cavity responses at carrier and sidebands
        F0 = self.reflection_coefficient(delta)
        Fp = self.reflection_coefficient(delta + self.eom.omega_m)
        Fm = self.reflection_coefficient(delta - self.eom.omega_m)

        # PDH baseband after mix + LPF
        Pc, Ps = self.laser.P_c, self.eom.P_s
        A = F0*np.conj(Fp) - np.conj(F0)*Fm
        Verr = 2*np.sqrt(Pc*Ps) * (np.real(A)*np.cos(demod_phase) + np.imag(A)*np.sin(demod_phase))
        return delta, Verr  

# PDH error signal 
rp = ReflectedPower()

# Time-averaged reflected power (no RF beat term)
delta = rp.omega_space - rp.laser.omega_c
F0 = rp.reflection_coefficient(delta)
Fp = rp.reflection_coefficient(delta + rp.eom.omega_m)
Fm = rp.reflection_coefficient(delta - rp.eom.omega_m)
P_avg = (rp.laser.P_c*np.abs(F0)**2
         + rp.eom.P_s*(np.abs(Fp)**2 + np.abs(Fm)**2))

delta, Verr = rp.error_signal(demod_phase=np.pi/2)

plt.figure(figsize=(9,6))
plt.subplot(2,1,1)
plt.plot(delta, P_avg)
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


