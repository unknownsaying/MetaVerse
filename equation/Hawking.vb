Imports System
Imports System.Collections.Generic
Imports System.Text

''' Black Hole equations:
''' Hawking Temperature: T = ħc³ / (8πGMk_B)
''' Bekenstein-Hawking Entropy: S = (k_B c³ A) / (4ħG)
''' Where A = 4πR_s² and R_s = 2GM/c²

Class HawkingBlackHoleThermodynamics
    ' Fundamental physical constants (SI units)
    Const Plank_constant As Double = 6.62607015e-34 ' J·s (h)
    Const reduced_Plank_constant As Double = 1.054571817e-34 ' ħ = h/(2π) J·s
    Const lightspeed As Double = 299792458.0 ' m/s
    Const gravity_constant As Double = 6.67430e-11 ' m³ kg⁻¹ s⁻²
    Const Boltzmann_constant As Double = 1.380649e-23 ' J/K
    Class SchwarzschildBlackHole
