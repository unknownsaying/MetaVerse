' MAXWELL'S EQUATIONS SUMMARY
' ===========================
' Differential Form (Microscopic)
Sub DifferentialForm()
    ' 1. Gauss's Law for Electricity: ∇·D = ρ_free
    ' 2. Gauss's Law for Magnetism: ∇·B = 0
    ' 3. Faraday's Law: ∇×E = -∂B/∂t
    ' 4. Ampère-Maxwell Law: ∇×H = J_free + ∂D/∂t
End Sub

' Integral Form (Macroscopic)
Sub IntegralForm()
    ' 1. ∮ D·dA = Q_free (Electric flux)
    ' 2. ∮ B·dA = 0 (Magnetic flux)
    ' 3. ∮ E·dl = -dΦ_B/dt (Faraday's Law)
    ' 4. ∮ H·dl = I_free + dΦ_D/dt (Ampère-Maxwell)
End Sub

' Constitutive Relations
Class ConstitutiveRelations
    Public D As String = "εE"  ' Permittivity: ε = ε₀ε_r
    Public B As String = "μH"  ' Permeability: μ = μ₀μ_r
    Public J As String = "σE"  ' Conductivity (Ohm's Law)
    
    Public ε0 As Double = 8.854187817e-12  ' F/m
    Public μ0 As Double = 1.2566370614e-6  ' N/A²
    Public c As Double = 299792458         ' m/s (c² = 1/(ε₀μ₀))
End Class

' Key Implications
Sub Implications()
    ' • Wave Equation: ∇²E = με ∂²E/∂t²
    ' • Continuity: ∇·J = -∂ρ/∂t
    ' • Lorentz Force: F = q(E + v×B)
    ' • Poynting Vector: S = E×H
    ' • Boundary Conditions: E_tangential, B_normal continuous
End Sub

' Potentials Formulation
Sub Potentials()
    ' E = -∇V - ∂A/∂t (V = scalar potential)
    ' B = ∇×A (A = vector potential)
    ' Lorenz Gauge: ∇·A + με ∂V/∂t = 0
    ' Wave equations for A and V
End Sub

' Special Cases
Sub SpecialCases()
    ' Electrostatics: ∇·D = ρ, ∇×E = 0
    ' Magnetostatics: ∇·B = 0, ∇×H = J
    ' Quasistatics: Time-varying, no radiation
    ' Vacuum: ρ=0, J=0 → EM waves at speed c
End Sub
' Main Summary
Sub Main()
    Console.WriteLine("MAXWELL'S EQUATIONS SUMMARY")
    Console.WriteLine("1. Electric charges create diverging E-fields (Gauss)")
    Console.WriteLine("2. No magnetic monopoles exist")
    Console.WriteLine("3. Changing B-fields induce circulating E-fields (Faraday)")
    Console.WriteLine("4. Currents & changing E-fields create circulating B-fields")
    Console.WriteLine("Key Unifications:")
    Console.WriteLine("- Unified electricity, magnetism, and optics")
    Console.WriteLine("- Predicted EM waves at speed of light")
    Console.WriteLine("- Showed light is electromagnetic wave")
    Console.WriteLine("- Basis for all classical electromagnetism")
End Sub