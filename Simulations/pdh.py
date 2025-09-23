from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt


class CONSTANTS:
    SPEED_OF_LIGHT = 3e8 # speed of light

@dataclass
class laser(CONSTANTS):
    l_0 : float # wavelength
    P_c : float # carrier power
    f_c : float = field(init=False) # carrier frequency
    omega_c : float = field(init=False)  # angular carrier frequency
    def __post_init__(self):
        self.omega_c = self.SPEED_OF_LIGHT/self.l_0
        self.f_c = self.omega_c/(2*np.pi)

@dataclass
class eom(CONSTANTS):
    f_m : float # modulation frequency
    omega_m : float # sideband angular frequency
    P_s : float # side band power
    
@dataclass 
class cavity(CONSTANTS):
    r : float  # reflection coefficient
    L : float  # length of the cavity
    fsr : float = field(init=False) # free spectral range
    def __post_init__(self):
        self.fsr = self.SPEED_OF_LIGHT/2*self.L

class ReflectionCoeff(cavity):
    def __init__(self,**kwargs):
        c =kwargs.get('cavity',cavity(r=0.95, L=0.02))
        super().__init__(r=c.r,L=c.L)

    # get reflection coeffient
    def reflection_coefficient(self,omega : np.linspace) -> float:
        n = self.r*(np.exp(1j*omega/self.fsr)-1)
        d = 1 - self.r*self.r*np.exp(1j*omega/self.fsr)
        return n/d

class ReflectedPower(ReflectionCoeff):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        
        self.laser = kwargs.get('laser',laser(l_0=1550e-9, P_c=5e-3))
        self.eom = kwargs.get('eom',eom(f_m=20e6, omega_m=2*np.pi*20e6, P_s=0.5e-3))
        self.omega_space = np.linspace(self.laser.omega_c-3*self.eom.omega_m, \
                                       self.laser.omega_c+3*self.eom.omega_m, \
                                       int(1e6))
        
    def reflected_power(self):
        t1 = np.real(self.laser.P_c*self.reflection_coefficient(self.omega_space - self.laser.omega_c)*np.conj(self.reflection_coefficient(self.omega_space - self.laser.omega_c)))
        
print(ReflectedPower().reflected_power())



