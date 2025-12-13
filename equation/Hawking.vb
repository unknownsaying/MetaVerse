Imports System
Imports System.Collections.Generic
Imports System.Text

''' <summary>
''' Implementation of Hawking's Black Hole thermodynamics equations:
''' 1. Hawking Temperature: T = ħc³ / (8πGMk_B)
''' 2. Bekenstein-Hawking Entropy: S = (k_B c³ A) / (4ħG)
''' Where A = 4πR_s² and R_s = 2GM/c²
''' </summary>
Public Class HawkingBlackHoleThermodynamics
    ' Fundamental physical constants (SI units)
    Public Const PLANCK_CONSTANT As Double = 6.62607015e-34 ' J·s (h)
    Public Const REDUCED_PLANCK_CONSTANT As Double = 1.054571817e-34 ' ħ = h/(2π) J·s
    Public Const SPEED_OF_LIGHT As Double = 299792458.0 ' m/s
    Public Const GRAVITATIONAL_CONSTANT As Double = 6.67430e-11 ' m³ kg⁻¹ s⁻²
    Public Const BOLTZMANN_CONSTANT As Double = 1.380649e-23 ' J/K
    Public Const PI As Double = Math.PI
    
    ' Quantum gravity scales
    Public Shared ReadOnly PLANCK_LENGTH As Double = Math.Sqrt(REDUCED_PLANCK_CONSTANT * GRAVITATIONAL_CONSTANT / Math.Pow(SPEED_OF_LIGHT, 3))
    Public Shared ReadOnly PLANCK_MASS As Double = Math.Sqrt(REDUCED_PLANCK_CONSTANT * SPEED_OF_LIGHT / GRAVITATIONAL_CONSTANT)
    Public Shared ReadOnly PLANCK_TEMPERATURE As Double = (Math.Pow(SPEED_OF_LIGHT, 5) * Math.Sqrt(REDUCED_PLANCK_CONSTANT / GRAVITATIONAL_CONSTANT)) / BOLTZMANN_CONSTANT
    Public Shared ReadOnly PLANCK_AREA As Double = 4 * PI * Math.Pow(PLANCK_LENGTH, 2)
    
    ''' <summary>
    ''' Represents a Schwarzschild black hole with mass M
    ''' </summary>
    Public Class SchwarzschildBlackHole
        Private _mass As Double ' kg
        Private _name As String
        
        Public Sub New(mass As Double, Optional name As String = "Black Hole")
            _mass = mass
            _name = name
        End Sub
        
        ''' <summary>
        ''' Schwarzschild radius: R_s = 2GM/c²
        ''' </summary>
        Public ReadOnly Property SchwarzschildRadius As Double
            Get
                Dim numerator As Double = 2 * GRAVITATIONAL_CONSTANT * _mass
                Dim denominator As Double = Math.Pow(SPEED_OF_LIGHT, 2)
                Return numerator / denominator
            End Get
        End Property
        
        ''' <summary>
        ''' Event horizon area: A = 4πR_s² = 16πG²M²/c⁴
        ''' </summary>
        Public ReadOnly Property HorizonArea As Double
            Get
                Dim rs As Double = SchwarzschildRadius
                Return 4 * PI * Math.Pow(rs, 2)
            End Get
        End Property
        
        ''' <summary>
        ''' Surface gravity at event horizon: κ = c⁴/(4GM)
        ''' </summary>
        Public ReadOnly Property SurfaceGravity As Double
            Get
                Dim numerator As Double = Math.Pow(SPEED_OF_LIGHT, 4)
                Dim denominator As Double = 4 * GRAVITATIONAL_CONSTANT * _mass
                Return numerator / denominator
            End Get
        End Property
        
        ''' <summary>
        ''' Hawking Temperature: T = ħc³/(8πGMk_B)
        ''' BREAKDOWN:
        '''   Numerator: ħc³
        '''   Denominator: 8πGMk_B
        ''' </summary>
        Public Function CalculateHawkingTemperature() As TemperatureComponents
            Dim components As New TemperatureComponents()
            
            ' Calculate numerator: ħc³
            components.Numerator = REDUCED_PLANCK_CONSTANT * Math.Pow(SPEED_OF_LIGHT, 3)
            components.NumeratorDescription = "ħc³ = " & FormatScientific(components.Numerator) & " J·m³/s³"
            components.NumeratorBreakdown = New Dictionary(Of String, Double) From {
                {"ħ (Reduced Planck Constant)", REDUCED_PLANCK_CONSTANT},
                {"c³ (Speed of Light Cubed)", Math.Pow(SPEED_OF_LIGHT, 3)}
            }
            
            ' Calculate denominator: 8πGMk_B
            components.Denominator = 8 * PI * GRAVITATIONAL_CONSTANT * _mass * BOLTZMANN_CONSTANT
            components.DenominatorDescription = "8πGMk_B = " & FormatScientific(components.Denominator) & " J·m²·kg/K"
            components.DenominatorBreakdown = New Dictionary(Of String, Double) From {
                {"8π", 8 * PI},
                {"G (Gravitational Constant)", GRAVITATIONAL_CONSTANT},
                {"M (Black Hole Mass)", _mass},
                {"k_B (Boltzmann Constant)", BOLTZMANN_CONSTANT}
            }
            
            ' Calculate final temperature
            components.Temperature = components.Numerator / components.Denominator
            components.TemperatureDescription = "T = ħc³/(8πGMk_B) = " & FormatScientific(components.Temperature) & " K"
            
            ' Alternative formula using surface gravity: T = ħκ/(2πk_Bc)
            Dim altTemperature As Double = (REDUCED_PLANCK_CONSTANT * SurfaceGravity) / (2 * PI * BOLTZMANN_CONSTANT * SPEED_OF_LIGHT)
            components.AlternativeForm = "T = ħκ/(2πk_Bc) = " & FormatScientific(altTemperature) & " K"
            
            Return components
        End Function
        
        ''' <summary>
        ''' Bekenstein-Hawking Entropy: S = (k_B c³ A)/(4ħG)
        ''' BREAKDOWN:
        '''   Numerator: k_B c³ A
        '''   Denominator: 4ħG
        ''' Equivalent form: S = (k_B A)/(4ℓ_P²) where ℓ_P is Planck length
        ''' </summary>
        Public Function CalculateBekensteinHawkingEntropy() As EntropyComponents
            Dim components As New EntropyComponents()
            Dim area As Double = HorizonArea
            
            ' Calculate numerator: k_B c³ A
            components.Numerator = BOLTZMANN_CONSTANT * Math.Pow(SPEED_OF_LIGHT, 3) * area
            components.NumeratorDescription = "k_B c³ A = " & FormatScientific(components.Numerator) & " J·m⁵/(K·s³)"
            components.NumeratorBreakdown = New Dictionary(Of String, Double) From {
                {"k_B (Boltzmann Constant)", BOLTZMANN_CONSTANT},
                {"c³ (Speed of Light Cubed)", Math.Pow(SPEED_OF_LIGHT, 3)},
                {"A (Horizon Area)", area}
            }
            
            ' Calculate denominator: 4ħG
            components.Denominator = 4 * REDUCED_PLANCK_CONSTANT * GRAVITATIONAL_CONSTANT
            components.DenominatorDescription = "4ħG = " & FormatScientific(components.Denominator) & " J·m³·kg/s"
            components.DenominatorBreakdown = New Dictionary(Of String, Double) From {
                {"4", 4.0},
                {"ħ (Reduced Planck Constant)", REDUCED_PLANCK_CONSTANT},
                {"G (Gravitational Constant)", GRAVITATIONAL_CONSTANT}
            }
            
            ' Calculate final entropy
            components.Entropy = components.Numerator / components.Denominator
            components.EntropyDescription = "S = (k_B c³ A)/(4ħG) = " & FormatScientific(components.Entropy) & " J/K"
            
            ' Alternative form using Planck length: S = (k_B A)/(4ℓ_P²)
            Dim altEntropy As Double = (BOLTZMANN_CONSTANT * area) / (4 * Math.Pow(PLANCK_LENGTH, 2))
            components.AlternativeForm = "S = (k_B A)/(4ℓ_P²) = " & FormatScientific(altEntropy) & " J/K"
            
            ' Calculate entropy in bits: S/k_B ln(2)
            components.EntropyInBits = components.Entropy / (BOLTZMANN_CONSTANT * Math.Log(2))
            components.EntropyBitsDescription = "S/(k_B ln2) = " & FormatNumber(components.EntropyInBits) & " bits"
            
            ' Calculate number of Planck areas
            components.PlanckAreas = area / PLANCK_AREA
            components.PlanckAreaDescription = $"A/A_P = {components.PlanckAreas:E6} Planck areas"
            
            Return components
        End Function
        
        ''' <summary>
        ''' Calculate black hole heat capacity: C = dE/dT = -8πGM²/ħc
        ''' </summary>
        Public ReadOnly Property HeatCapacity As Double
            Get
                Dim numerator As Double = -8 * PI * GRAVITATIONAL_CONSTANT * Math.Pow(_mass, 2)
                Dim denominator As Double = REDUCED_PLANCK_CONSTANT * SPEED_OF_LIGHT
                Return numerator / denominator
            End Get
        End Property
        
        ''' <summary>
        ''' Calculate black hole lifetime via Hawking radiation
        ''' τ ~ (5120πG²M³)/(ħc⁴)
        ''' </summary>
        Public ReadOnly Property HawkingLifetime As Double
            Get
                Dim numerator As Double = 5120 * PI * Math.Pow(GRAVITATIONAL_CONSTANT, 2) * Math.Pow(_mass, 3)
                Dim denominator As Double = REDUCED_PLANCK_CONSTANT * Math.Pow(SPEED_OF_LIGHT, 4)
                Return numerator / denominator
            End Get
        End Property
        
        ''' <summary>
        ''' Black hole mass-energy equivalence: E = Mc²
        ''' </summary>
        Public ReadOnly Property TotalEnergy As Double
            Get
                Return _mass * Math.Pow(SPEED_OF_LIGHT, 2)
            End Get
        End Property
        
        Public Property Mass As Double
            Get
                Return _mass
            End Get
            Set(value As Double)
                _mass = value
            End Set
        End Property
        
        Public ReadOnly Property Name As String
            Get
                Return _name
            End Get
        End Property
        
        Public Overrides Function ToString() As String
            Return $"{_name}: M = {FormatScientific(_mass)} kg, R_s = {FormatScientific(SchwarzschildRadius)} m"
        End Function
    End Class
    
    ''' <summary>
    ''' Detailed breakdown of Hawking Temperature calculation
    ''' </summary>
    Public Class TemperatureComponents
        Public Property Temperature As Double
        Public Property Numerator As Double
        Public Property Denominator As Double
        Public Property NumeratorDescription As String
        Public Property DenominatorDescription As String
        Public Property TemperatureDescription As String
        Public Property AlternativeForm As String
        Public Property NumeratorBreakdown As Dictionary(Of String, Double)
        Public Property DenominatorBreakdown As Dictionary(Of String, Double)
        
        Public Sub New()
            NumeratorBreakdown = New Dictionary(Of String, Double)()
            DenominatorBreakdown = New Dictionary(Of String, Double)()
        End Sub
        
        Public Function GetDetailedReport() As String
            Dim sb As New StringBuilder()
            sb.AppendLine("HAWKING TEMPERATURE BREAKDOWN:")
            sb.AppendLine("Formula: T = ħc³ / (8πGMk_B)")
            sb.AppendLine()
            
            sb.AppendLine("NUMERATOR (ħc³):")
            For Each kvp In NumeratorBreakdown
                sb.AppendLine($"  {kvp.Key}: {FormatScientific(kvp.Value)}")
            Next
            sb.AppendLine($"  Total Numerator: {FormatScientific(Numerator)}")
            sb.AppendLine($"  Description: {NumeratorDescription}")
            sb.AppendLine()
            
            sb.AppendLine("DENOMINATOR (8πGMk_B):")
            For Each kvp In DenominatorBreakdown
                sb.AppendLine($"  {kvp.Key}: {FormatScientific(kvp.Value)}")
            Next
            sb.AppendLine($"  Total Denominator: {FormatScientific(Denominator)}")
            sb.AppendLine($"  Description: {DenominatorDescription}")
            sb.AppendLine()
            
            sb.AppendLine("FINAL CALCULATION:")
            sb.AppendLine($"  T = {FormatScientific(Numerator)} / {FormatScientific(Denominator)}")
            sb.AppendLine($"  T = {FormatScientific(Temperature)} K")
            sb.AppendLine()
            sb.AppendLine("ALTERNATIVE FORM:")
            sb.AppendLine($"  {AlternativeForm}")
            
            Return sb.ToString()
        End Function
    End Class
    
    ''' <summary>
    ''' Detailed breakdown of Bekenstein-Hawking Entropy calculation
    ''' </summary>
    Public Class EntropyComponents
        Public Property Entropy As Double
        Public Property Numerator As Double
        Public Property Denominator As Double
        Public Property NumeratorDescription As String
        Public Property DenominatorDescription As String
        Public Property EntropyDescription As String
        Public Property AlternativeForm As String
        Public Property EntropyInBits As Double
        Public Property EntropyBitsDescription As String
        Public Property PlanckAreas As Double
        Public Property PlanckAreaDescription As String
        Public Property NumeratorBreakdown As Dictionary(Of String, Double)
        Public Property DenominatorBreakdown As Dictionary(Of String, Double)
        
        Public Sub New()
            NumeratorBreakdown = New Dictionary(Of String, Double)()
            DenominatorBreakdown = New Dictionary(Of String, Double)()
        End Sub
        
        Public Function GetDetailedReport() As String
            Dim sb As New StringBuilder()
            sb.AppendLine("BEKENSTEIN-HAWKING ENTROPY BREAKDOWN:")
            sb.AppendLine("Formula: S = (k_B c³ A) / (4ħG)")
            sb.AppendLine()
            
            sb.AppendLine("NUMERATOR (k_B c³ A):")
            For Each kvp In NumeratorBreakdown
                sb.AppendLine($"  {kvp.Key}: {FormatScientific(kvp.Value)}")
            Next
            sb.AppendLine($"  Total Numerator: {FormatScientific(Numerator)}")
            sb.AppendLine($"  Description: {NumeratorDescription}")
            sb.AppendLine()
            
            sb.AppendLine("DENOMINATOR (4ħG):")
            For Each kvp In DenominatorBreakdown
                sb.AppendLine($"  {kvp.Key}: {FormatScientific(kvp.Value)}")
            Next
            sb.AppendLine($"  Total Denominator: {FormatScientific(Denominator)}")
            sb.AppendLine($"  Description: {DenominatorDescription}")
            sb.AppendLine()
            
            sb.AppendLine("FINAL CALCULATION:")
            sb.AppendLine($"  S = {FormatScientific(Numerator)} / {FormatScientific(Denominator)}")
            sb.AppendLine($"  S = {FormatScientific(Entropy)} J/K")
            sb.AppendLine($"  {EntropyDescription}")
            sb.AppendLine()
            
            sb.AppendLine("QUANTUM INFORMATION INTERPRETATION:")
            sb.AppendLine($"  Entropy in bits: {EntropyBitsDescription}")
            sb.AppendLine($"  Planck areas: {PlanckAreaDescription}")
            sb.AppendLine($"  Information content: ~{FormatNumber(EntropyInBits)} qubits")
            sb.AppendLine()
            
            sb.AppendLine("ALTERNATIVE FORM (Using Planck Length):")
            sb.AppendLine($"  {AlternativeForm}")
            
            Return sb.ToString()
        End Function
    End Class
    
    ''' <summary>
    ''' Demonstrates the holographic principle: maximum entropy in a region
    ''' is proportional to its surface area, not volume
    ''' </summary>
    Public Class HolographicPrinciple
        ''' <summary>
        ''' Calculate maximum entropy for a spherical region of radius R
        ''' S_max = (k_B c³ A)/(4ħG) = (πk_B R²)/(ħG/c³)
        ''' </summary>
        Public Shared Function CalculateMaximumEntropy(radius As Double) As Double
            Dim area As Double = 4 * PI * Math.Pow(radius, 2)
            Return (BOLTZMANN_CONSTANT * area) / (4 * Math.Pow(PLANCK_LENGTH, 2))
        End Function
        
        ''' <summary>
        ''' Calculate Bekenstein bound: maximum information in a region
        ''' I_max = (2πER)/(ħc ln2)
        ''' </summary>
        Public Shared Function CalculateBekensteinBound(energy As Double, radius As Double) As Double
            Dim numerator As Double = 2 * PI * energy * radius
            Dim denominator As Double = REDUCED_PLANCK_CONSTANT * SPEED_OF_LIGHT * Math.Log(2)
            Return numerator / denominator
        End Function
    End Class
    
    ''' <summary>
    ''' Helper function to format numbers in scientific notation
    ''' </summary>
    Private Shared Function FormatScientific(value As Double) As String
        Return String.Format("{0:E6}", value)
    End Function
    
    ''' <summary>
    ''' Helper function to format large numbers
    ''' </summary>
    Private Shared Function FormatNumber(value As Double) As String
        If value >= 1e100 Then
            Return String.Format("{0:E6}", value)
        ElseIf value >= 1e6 Then
            Return String.Format("{0:E6}", value)
        Else
            Return String.Format("{0:F2}", value)
        End If
    End Function
    
    ''' <summary>
    ''' Calculate Planck scales for reference
    ''' </summary>
    Public Shared Function GetPlanckScalesReport() As String
        Dim sb As New StringBuilder()
        sb.AppendLine("PLANCK SCALES (Quantum Gravity):")
        sb.AppendLine($"  Planck Length (ℓ_P): {FormatScientific(PLANCK_LENGTH)} m")
        sb.AppendLine($"  Planck Mass (m_P): {FormatScientific(PLANCK_MASS)} kg")
        sb.AppendLine($"  Planck Temperature (T_P): {FormatScientific(PLANCK_TEMPERATURE)} K")
        sb.AppendLine($"  Planck Area (A_P): {FormatScientific(PLANCK_AREA)} m²")
        Return sb.ToString()
    End Function
End Class

''' <summary>
''' Main demonstration of Hawking Black Hole Thermodynamics
''' </summary>
Module HawkingThermodynamicsDemo
    Sub Main()
        Console.WriteLine("HAWKING BLACK HOLE THERMODYNAMICS IN VB.NET")
        Console.WriteLine("=".PadRight(70, "="))
        Console.WriteLine()
        
        ' Display Planck scales
        Console.WriteLine(HawkingBlackHoleThermodynamics.GetPlanckScalesReport())
        
        ' Example 1: Solar mass black hole
        Console.WriteLine("EXAMPLE 1: SOLAR MASS BLACK HOLE")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Dim solarMass As Double = 1.989e30 ' kg
        Dim blackHole1 As New HawkingBlackHoleThermodynamics.SchwarzschildBlackHole(solarMass, "Solar Mass BH")
        
        Console.WriteLine(blackHole1.ToString())
        Console.WriteLine($"  Horizon Area: {FormatScientific(blackHole1.HorizonArea)} m²")
        Console.WriteLine($"  Surface Gravity: {FormatScientific(blackHole1.SurfaceGravity)} m/s²")
        Console.WriteLine($"  Total Energy: {FormatScientific(blackHole1.TotalEnergy)} J")
        
        Dim temp1 = blackHole1.CalculateHawkingTemperature()
        Console.WriteLine()
        Console.WriteLine(temp1.GetDetailedReport())
        
        Dim entropy1 = blackHole1.CalculateBekensteinHawkingEntropy()
        Console.WriteLine(entropy1.GetDetailedReport())
        
        Console.WriteLine($"Heat Capacity: {FormatScientific(blackHole1.HeatCapacity)} J/K")
        Console.WriteLine($"Hawking Lifetime: {FormatScientific(blackHole1.HawkingLifetime)} s")
        Console.WriteLine($"  ≈ {FormatNumber(blackHole1.HawkingLifetime / (365.25 * 24 * 3600 * 1e9))} billion years")
        
        ' Example 2: Planck mass black hole (quantum black hole)
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 2: PLANCK MASS BLACK HOLE (QUANTUM)")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Dim planckMassBlackHole As New HawkingBlackHoleThermodynamics.SchwarzschildBlackHole(
            HawkingBlackHoleThermodynamics.PLANCK_MASS, "Planck Mass BH")
        
        Console.WriteLine(planckMassBlackHole.ToString())
        
        Dim temp2 = planckMassBlackHole.CalculateHawkingTemperature()
        Dim entropy2 = planckMassBlackHole.CalculateBekensteinHawkingEntropy()
        
        Console.WriteLine($"Hawking Temperature: {FormatScientific(temp2.Temperature)} K")
        Console.WriteLine($"Entropy: {FormatScientific(entropy2.Entropy)} J/K")
        Console.WriteLine($"Entropy in bits: {FormatNumber(entropy2.EntropyInBits)}")
        Console.WriteLine($"Heat Capacity: {FormatScientific(planckMassBlackHole.HeatCapacity)} J/K")
        Console.WriteLine($"Hawking Lifetime: {FormatScientific(planckMassBlackHole.HawkingLifetime)} s")
        
        ' Example 3: Mini black hole (1 billion tons)
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 3: MINI BLACK HOLE (1 BILLION TONS)")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Dim miniBlackHoleMass As Double = 1e12 ' 1 billion tons = 1e12 kg
        Dim miniBlackHole As New HawkingBlackHoleThermodynamics.SchwarzschildBlackHole(miniBlackHoleMass, "Mini BH")
        
        Console.WriteLine(miniBlackHole.ToString())
        
        Dim temp3 = miniBlackHole.CalculateHawkingTemperature()
        Dim entropy3 = miniBlackHole.CalculateBekensteinHawkingEntropy()
        
        Console.WriteLine($"Hawking Temperature: {FormatScientific(temp3.Temperature)} K")
        Console.WriteLine($"Entropy: {FormatScientific(entropy3.Entropy)} J/K")
        Console.WriteLine($"Entropy in bits: {FormatNumber(entropy3.EntropyInBits)}")
        Console.WriteLine($"Hawking Lifetime: {FormatScientific(miniBlackHole.HawkingLifetime)} s")
        Console.WriteLine($"  ≈ {FormatNumber(miniBlackHole.HawkingLifetime)} seconds")
        
        ' Example 4: Supermassive black hole (Sagittarius A*)
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 4: SUPERMASSIVE BLACK HOLE (SAGITTARIUS A*)")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Dim sagittariusAMass As Double = 4.154e36 ' kg (≈4.154 million solar masses)
        Dim sagittariusA As New HawkingBlackHoleThermodynamics.SchwarzschildBlackHole(sagittariusAMass, "Sagittarius A*")
        
        Console.WriteLine(sagittariusA.ToString())
        
        Dim temp4 = sagittariusA.CalculateHawkingTemperature()
        Dim entropy4 = sagittariusA.CalculateBekensteinHawkingEntropy()
        
        Console.WriteLine($"Hawking Temperature: {FormatScientific(temp4.Temperature)} K")
        Console.WriteLine($"Entropy: {FormatScientific(entropy4.Entropy)} J/K")
        Console.WriteLine($"Entropy in bits: {FormatNumber(entropy4.EntropyInBits)}")
        Console.WriteLine($"  (Compare: Observable universe entropy ~10^{120} bits)")
        
        ' Example 5: Holographic principle demonstration
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 5: HOLOGRAPHIC PRINCIPLE")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Dim radius As Double = 1e26 ' Approximate observable universe radius (m)
        Dim maxEntropy = HawkingBlackHoleThermodynamics.HolographicPrinciple.CalculateMaximumEntropy(radius)
        Dim universeArea As Double = 4 * Math.PI * Math.Pow(radius, 2)
        
        Console.WriteLine($"Observable Universe Radius: {FormatScientific(radius)} m")
        Console.WriteLine($"Surface Area: {FormatScientific(universeArea)} m²")
        Console.WriteLine($"Maximum Entropy (Holographic Bound): {FormatScientific(maxEntropy)} J/K")
        Console.WriteLine($"Maximum Information: ~{FormatNumber(maxEntropy / (HawkingBlackHoleThermodynamics.BOLTZMANN_CONSTANT * Math.Log(2)))} bits")
        
        ' Example 6: Black hole thermodynamics relationships
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 6: THERMODYNAMIC RELATIONSHIPS")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Console.WriteLine("First Law of Black Hole Thermodynamics:")
        Console.WriteLine("  dM = T dS + Ω dJ + Φ dQ")
        Console.WriteLine("  Where: M = mass, T = temperature, S = entropy,")
        Console.WriteLine("         Ω = angular velocity, J = angular momentum,")
        Console.WriteLine("         Φ = electric potential, Q = electric charge")
        Console.WriteLine()
        
        Console.WriteLine("Second Law (Generalized):")
        Console.WriteLine("  δS_total = δS_BH + δS_matter ≥ 0")
        Console.WriteLine("  Black hole area never decreases in classical GR")
        Console.WriteLine()
        
        Console.WriteLine("Third Law (Nernst):")
        Console.WriteLine("  T → 0 as surface gravity κ → 0")
        Console.WriteLine("  Extremal black holes have zero temperature")
        
        ' Example 7: Entropy vs Temperature relationship
        Console.WriteLine()
        Console.WriteLine("EXAMPLE 7: ENTROPY-TEMPERATURE RELATIONSHIP")
        Console.WriteLine("-".PadRight(60, "-"))
        
        Console.WriteLine("For Schwarzschild black holes:")
        Console.WriteLine("  S ∝ M²  and  T ∝ 1/M")
        Console.WriteLine("  Therefore: S ∝ 1/T²")
        Console.WriteLine("  This shows negative heat capacity: C = dM/dT < 0")
        
        ' Summary of key equations
        Console.WriteLine()
        Console.WriteLine("SUMMARY OF KEY EQUATIONS:")
        Console.WriteLine("=".PadRight(60, "="))
        Console.WriteLine()
        Console.WriteLine("1. Schwarzschild Radius:")
        Console.WriteLine("   R_s = 2GM/c²")
        Console.WriteLine()
        Console.WriteLine("2. Hawking Temperature:")
        Console.WriteLine("   T = ħc³/(8πGMk_B)")
        Console.WriteLine("     = ħκ/(2πk_Bc)")
        Console.WriteLine()
        Console.WriteLine("3. Bekenstein-Hawking Entropy:")
        Console.WriteLine("   S = (k_B c³ A)/(4ħG)")
        Console.WriteLine("     = (k_B A)/(4ℓ_P²)")
        Console.WriteLine("     ≈ 1.4 × 10⁶⁹ (M/M☉)² J/K")
        Console.WriteLine()
        Console.WriteLine("4. Thermodynamic Relations:")
        Console.WriteLine("   dM = T dS  (for non-rotating, uncharged)")
        Console.WriteLine("   C = dM/dT = -8πGM²/ħc")
        Console.WriteLine("   F = M - TS = M/2  (Helmholtz free energy)")
        
        Console.WriteLine()
        Console.WriteLine("Press any key to exit...")
        Console.ReadKey()
    End Sub
    
    Private Function FormatScientific(value As Double) As String
        Return String.Format("{0:E6}", value)
    End Function
    
    Private Function FormatNumber(value As Double) As String
        Return String.Format("{0:E6}", value)
    End Function
End Module