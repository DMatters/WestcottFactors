from .user import *

class Kinematics(UserSpectrum):
    __doc__="""Class to handle quantities related to neutron-beam kinematics."""

    # Constants
    eV = 1.602189e-19 #J
    kB = 1.38066e-23 #Boltzmann's constant, J/K
    m_n = 1.00866501 *1.660566e-27 #neutron mass, kg
    kB_eVK = 8.6117343 #Boltzmann's constant, eV/K
    h = 6.626183e-34 #Planck's constant, J*s
    h_eVs = 4.13567e-15 #Planck's constant, eV*s
   
    # Thermal quantities
    v_0 = 2200 #thermal-neutron velocity, m/s
    E_0 = 1/2 * m_n * v_0**2 / eV #eV
    T_0 = 293 # K
    lambda_0 = 1.8 # A

    def __init__(self):
        pass

    def display_thermal_properties(self) -> None:
        """Prints to screen thermal neutron properties:
        -Velocity
        -Temperature
        -Wavelength
        -Energy"""
        print(f"Velocity = {Kinematics.v_0} m/s")
        print(f"Temperature = {Kinematics.T_0} K")
        print(f"Wavelength = {Kinematics.lambda_0} A")
        print(f"Energy = {1000*Kinematics.E_0} meV")

    def display_constants(self) -> None:
        """Prints to screen various constants useful for dealing with neutron 
        beam kinematics and unit conversions."""
        print(f"1 eV = {Kinematics.eV} J")
        print(f"Neutron mass = {Kinematics.m_n} kg")
        print(f"Boltzmann constant = {Kinematics.kB} J/K")
        print(f"Boltzmann constant = {Kinematics.kB_eVK} eV/K")
        print(f"Planck constant = {Kinematics.h} J*s")
        print(f"Planck constant = {Kinematics.h_eVs} eV*s")
        
    def vel(self,E):
        """Convert neutron energy (eV) to velocity (m/s)"""
        self.E = E
        
        E_joules = self.E*Kinematics.eV
        return np.sqrt(2*E_joules/Kinematics.m_n)

    def phi_Maxwellian(self, T, v_array):
        """Maxwellian velocity distribution at a given temperature T (K)"""
        self.T = T
        self.v_array = np.array(v_array)
        
        phi = []
        vt = np.sqrt(2*Kinematics.kB*self.T/Kinematics.m_n)
        for v in self.v_array:
            phi.append(2 * np.exp(-v**2/vt**2) * v**3/vt**4)
        return np.array(phi)


class Irregularity(Kinematics):
    __doc__="""Class to handle irregularity methods in the approximation of the 
    Westcott g-factor calculation."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def del_0(self, v, E_resonance, Gamma):
        """Lorentzian lineshape for the irregularity function, per Molnar 
        Eqs. 1-3."""
        self.v = v
        self.E_resonance = E_resonance
        self.Gamma = Gamma
        
        K = Kinematics()
        E = 0.5 * K.m_n * self.v**2 / K.eV
        return ((self.E_resonance - K.E_0)**2 + (self.Gamma**2)/4) / ((self.E_resonance - E)**2 + (self.Gamma**2)/4)

    def p(self, vn, T=None):
        """Neutron density function (Molnar Ch. 1, Table 1)."""
        self.vn = vn
        self.T = T

        K = Kinematics()
        phi = None
        if T is None:
            phi = K.phi_Maxwellian(K.T_0, np.logspace(0,5,100000))
        else:
            phi = K.phi_Maxwellian(self.T, self.vn)

        p_array = []
        for i in range(len(self.vn)):
            vt = np.sqrt(2 * K.kB * self.T / K.m_n)
            p_array.append(2 * vt * phi[i] / (np.sqrt(np.pi) * self.vn[i]))
        N = trapezoid(np.array(p_array), self.vn)  # Normalization factor, to ensure integral of p(T,v) integrates to unity (Molnar p. 12)
        return np.array(p_array)/N

    def gw_irregularity(self, E_resonance, Gamma, T=None, vn=np.logspace(0,5,100000)):
        """Evaluate g-factor using irregularity function method described by 
        Molnar et al. (Eqs. 1-5)."""
        self.E_resonance = E_resonance
        self.Gamma = Gamma
        self.T = T
        self.vn = vn

        d = []
        for i in range(len(self.vn)):
            d.append(Irregularity.del_0(self, self.vn[i], self.E_resonance, self.Gamma))
        d = np.array(d)
        return trapezoid(d * Irregularity.p(self, self.vn, self.T), self.vn)
    
class gFactors(Irregularity):
    __doc__="""Class to handle numerical integration of complete cross-section 
    spectrum for the determination of Westcott g-factors."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def gw_Maxwellian(self, T, E, sigma, vn=np.logspace(0,5,100000)):
        """Westcott g-factor according to assumed theoretical Maxwellian 
        distribution at a given neutron temperature."""
        self.T = T
        self.E = E
        self.sigma = sigma

        K = Kinematics()
        v_sigma = K.vel(self.E)  #convert neutron energy to velocity
        dndv = K.phi_Maxwellian(self.T, vn)
        sigma0 = np.interp(K.v_0, v_sigma, self.sigma)  #thermal cross section, barns
        sigma_interp = np.interp(vn, v_sigma, self.sigma)

        return 1/(sigma0 * K.v_0) * trapezoid(dndv * vn * sigma_interp, vn) / trapezoid(dndv, vn)

    def gw_arbitrary(self, E_spectrum, dndE_spectrum, E_endf, sigma_endf, vn=np.logspace(0,5,100000)):
        """Integrate to evaluate Westcott g-factor for an arbitrary neutron 
        flux distribution."""
        self.E_spectrum = E_spectrum
        self.dndE_spectrum = dndE_spectrum
        self.E_endf = E_endf
        self.sigma_endf = sigma_endf

        K = Kinematics()
        E_n = 0.5*K.m_n * vn**2/K.eV  #energy space, eV
        dndE_interp = np.interp(E_n, self.E_spectrum, self.dndE_spectrum, left=0, right=0)
        dndv_interp = np.sqrt(2 * K.m_n * E_n) * dndE_interp
        v_sigma = K.vel(self.E_endf)
        sigma_interp = np.interp(vn, v_sigma, self.sigma_endf)
        sigma0 = np.interp(K.v_0, v_sigma, self.sigma_endf)  #thermal cross section, barns
    
        return 1/(sigma0 * K.v_0) * trapezoid(dndv_interp * vn * sigma_interp, vn) / trapezoid(dndv_interp, vn)

    
