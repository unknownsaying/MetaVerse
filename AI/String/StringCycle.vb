Public Class GuitarString
    Public Property Note As String
    Public Property Frequency As Double ' in Hz
    Public Property StringColor As Color
    Public Property YPosition As Integer
    Public Property Length As Integer
    Public Property Tension As Double
    Public Property MassDensity As Double ' μ in physics equations
    
    ' String vibration state
    Public Property IsPlucked As Boolean = False
    Public Property PluckTime As Double = 0
    
    Public Overrides Function ToString() As String
        Return $"{Note} String: f={Frequency}Hz, T={Tension}, μ={MassDensity}"
    End Function
End Class

Public Module StringPhysics
    ' Calculate wave speed on string
    Public Function WaveSpeed(tension As Double, massDensity As Double) As Double
        ' v = √(T/μ)
        Return Math.Sqrt(tension / massDensity)
    End Function

    ' Calculate harmonic frequencies
    Public Function HarmonicFrequencies(fundamental As Double, 
                                        harmonicCount As Integer) As Double()
        Dim harmonics(harmonicCount - 1) As Double
        For i As Integer = 0 To harmonicCount - 1
            harmonics(i) = fundamental * (i + 1)
        Next
        Return harmonics
    End Function

    ' Superstring vibration modes (inspired by string theory)
    Public Function SuperstringMode(n As Integer, 
                                    dimension As Integer, 
                                    tension As Double) As Double
        ' Simplified version of string theory vibration energy
        Return Math.Sqrt(n * tension * (1 + dimension / 10.0))
    End Function

    ' Calculate standing wave pattern
    Public Function StandingWave(x As Double, 
                                 t As Double, 
                                 frequency As Double, 
                                 length As Double, 
                                 mode As Integer) As Double
        Dim k As Double = mode * Math.PI / length
        Dim omega As Double = 2 * Math.PI * frequency
        Return Math.Sin(k * x) * Math.Cos(omega * t)
    End Function
End Module

Imports System
Imports System.Windows.Forms

Module Program
    <STAThread>
    Sub Main()
        Application.EnableVisualStyles()
        Application.SetCompatibleTextRenderingDefault(False)
        Application.Run(New MainForm())
    End Sub
End Module