Class GuitarString
    Property Note As String
    Property Frequency As Double ' in Hz
    Property StringColor As Color
    Property YPosition As Integer
    Property Length As Integer
    Property Tension As Double
    Property MassDensity As Double ' μ in physics equations
    
    ' String vibration state
    Property IsPlucked As Boolean = False
    Property PluckTime As Double = 0
    
    Overrides Function ToString() As String
        Return $"{Note} String: f={Frequency}Hz, T={Tension}, μ={MassDensity}"
    End Function
End Class

Module StringPhysics
    ' Calculate wave speed on string
    Function WaveSpeed(tension As Double, massDensity As Double) As Double
        ' v = √(T/μ)
        Return Math.Sqrt(tension / massDensity)
    End Function

    ' Calculate harmonic frequencies
    Function HarmonicFrequencies(fundamental As Double, 
                                        harmonicCount As Integer) As Double()
        Dim harmonics(harmonicCount - 1) As Double
        For i As Integer = 0 To harmonicCount - 1
            harmonics(i) = fundamental * (i + 1)
        Next
        Return harmonics
    End Function

    ' Superstring vibration modes (inspired by string theory)
    Function SuperstringMode(n As Integer, 
                                    dimension As Integer, 
                                    tension As Double) As Double
        ' Simplified version of string theory vibration energy
        Return Math.Sqrt(n * tension * (1 + dimension / 10.0))
    End Function

    ' Calculate standing wave pattern
    Function StandingWave(x As Double, 
                                 t As Double, 
                                 frequency As Double, 
                                 length As Double, 
                                 mode As Integer) As Double
        Dim k As Double = mode * Math.PI / length
        Dim omega As Double = 2 * Math.PI * frequency
        Return Math.Sin(k * x) * Math.Cos(omega * t)
    End Function
End Module

Module Program
    <STAThread>
    Sub Main()
        Application.EnableVisualStyles()
        Application.SetCompatibleTextRenderingDefault(False)
        Application.Run(New MainForm())
    End Sub

End Module
