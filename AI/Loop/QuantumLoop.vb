Imports System.Drawing.Drawing2D
Imports System.Numerics

Class QuantumPhysicsForm
    Declare WithEvents simulationTimer As New Timer()
    Declare quantumParticles As New List(Of QuantumParticle)
    Declare entangledPairs As New List(Of EntangledPair)
    Declare superpositionParticles As New List(Of SuperpositionState)
    Declare waveFunction As New WaveFunction()
    Declare observer As New QuantumObserver()
    Declare measurementHistory As New List(Of MeasurementEvent)
    Declare rnd As New Random()
    Declare simulationTime As Double = 0
    Declare isCollapsing As Boolean = False

    ' Quantum Constants
    Declare Const PLANCK_CONSTANT As Double = 1.0545718E-34
    Declare Const SPEED_OF_LIGHT As Double = 299792458
    Declare Const BOLTZMANN_CONSTANT As Double = 1.380649E-23

    Class QuantumParticle
        Property Position As Vector2
        Property Momentum As Vector2
        Property Spin As Double ' ±1/2
        Property Charge As Double ' ±1
        Property Mass As Double
        Property WaveFunction As Complex ' Ψ(x,t)
        Property ProbabilityAmplitude As Double
        Property QuantumState As Integer ' 0 or 1
        Property IsEntangled As Boolean
        Property EntanglementID As Guid
        Property DecoherenceTime As Double
        Property Color As Color
        Property Size As Single
        Property EnergyLevel As Integer

        Sub New(pos As Vector2, mom As Vector2, spin As Double, chg As Double)
            Position = pos
            Momentum = mom
            Spin = spin
            Charge = chg
            Mass = 9.10938356E-31 ' Electron mass
            WaveFunction = Complex.One
            ProbabilityAmplitude = 1.0
            QuantumState = If(rnd.NextDouble() > 0.5, 1, 0)
            DecoherenceTime = rnd.NextDouble() * 10
            Color = Color.FromArgb(200, 100, 200, 255)
            Size = 6.0F
            EnergyLevel = 1
        End Sub

        Function CalculateDeBroglieWavelength() As Double
            If Momentum.Length() = 0 Then Return Double.MaxValue
            Return (2 * Math.PI * PLANCK_CONSTANT) / (Mass * Momentum.Length())
        End Function

        Sub ApplyUncertaintyPrinciple(uncertainty As Double)
            ' Δx * Δp ≥ ħ/2
            Dim deltaX As Double = rnd.NextDouble() * uncertainty
            Dim deltaP As Double = PLANCK_CONSTANT / (2 * deltaX)
            
            Position += New Vector2(CSng(deltaX * (rnd.NextDouble() - 0.5)), 
                                   CSng(deltaX * (rnd.NextDouble() - 0.5)))
            Momentum += New Vector2(CSng(deltaP * (rnd.NextDouble() - 0.5)), 
                                   CSng(deltaP * (rnd.NextDouble() - 0.5)))
        End Sub
    End Class

    Class SuperpositionState
        Property BaseParticle As QuantumParticle
        Property PossibleStates As List(Of Vector2)
        Property ProbabilityDistribution As List(Of Double)
        Property PhaseAngles As List(Of Double)
        Property IsCollapsed As Boolean = False

        Sub New(particle As QuantumParticle)
            BaseParticle = particle
            PossibleStates = New List(Of Vector2)
            ProbabilityDistribution = New List(Of Double)
            PhaseAngles = New List(Of Double)
            
            ' Create superposition of 3 possible states
            For i As Integer = 0 To 2
                PossibleStates.Add(particle.Position + 
                    New Vector2(CSng((rnd.NextDouble() - 0.5) * 100),
                               CSng((rnd.NextDouble() - 0.5) * 100)))
                
                Dim prob As Double = 1.0 / 3.0
                ProbabilityDistribution.Add(prob)
                PhaseAngles.Add(rnd.NextDouble() * Math.PI * 2)
            Next
        End Sub

        Function CalculateWaveFunction(x As Vector2, t As Double) As Complex
            Dim result As Complex = Complex.Zero
            
            For i As Integer = 0 To PossibleStates.Count - 1
                Dim distance As Double = Vector2.Distance(x, PossibleStates(i))
                Dim phase As Double = PhaseAngles(i) - EnergyAtState(i) * t / PLANCK_CONSTANT
                
                ' ψ(x,t) = A * e^(i(kx - ωt))
                Dim amplitude As Double = ProbabilityDistribution(i)
                Dim wave As New Complex(
                    amplitude * Math.Cos(phase),
                    amplitude * Math.Sin(phase))
                    
                ' Gaussian wave packet
                Dim gaussian As Double = Math.Exp(-distance * distance / 50.0)
                result += wave * gaussian
            Next
            
            Return result
        End Function

        Declare Function EnergyAtState(stateIndex As Integer) As Double
            ' E = ħω
            Return (stateIndex + 1) * 1.602E-19 ' Simple energy levels in eV
        End Function
    End Class

    Class EntangledPair
        Property ParticleA As QuantumParticle
        Property ParticleB As QuantumParticle
        Property CorrelationType As String ' "spin", "position", "momentum"
        Property BellState As Integer ' 0: Φ⁺, 1: Φ⁻, 2: Ψ⁺, 3: Ψ⁻
        Property EntanglementStrength As Double ' 0 to 1

        Sub New(p1 As QuantumParticle, p2 As QuantumParticle, corrType As String)
            ParticleA = p1
            ParticleB = p2
            CorrelationType = corrType
            BellState = rnd.Next(0, 4)
            EntanglementStrength = 0.95 + rnd.NextDouble() * 0.05
            
            ' Make particles entangled
            p1.IsEntangled = True
            p2.IsEntangled = True
            p1.EntanglementID = Guid.NewGuid()
            p2.EntanglementID = p1.EntanglementID
            
            ' Initialize with perfect correlation
            If CorrelationType = "spin" Then
                If BellState = 0 Or BellState = 1 Then
                    p2.Spin = p1.Spin ' Same spin
                Else
                    p2.Spin = -p1.Spin ' Opposite spin
                End If
            End If
        End Sub

        Sub UpdateCorrelation()
            ' Maintain entanglement through interaction
            If CorrelationType = "spin" Then
                If BellState = 0 Or BellState = 1 Then
                    ParticleB.Spin = ParticleA.Spin * CSng(EntanglementStrength)
                Else
                    ParticleB.Spin = -ParticleA.Spin * CSng(EntanglementStrength)
                End If
            ElseIf CorrelationType = "position" Then
                ' Position entanglement (EPR paradox)
                Dim delta As Vector2 = ParticleA.Position - ParticleB.Position
                If delta.Length() > 100 Then
                    ParticleB.Position = ParticleA.Position + 
                        New Vector2(CSng((rnd.NextDouble() - 0.5) * 10),
                                   CSng((rnd.NextDouble() - 0.5) * 10))
                End If
            End If
        End Sub
    End Class

    Class WaveFunction
        Property PotentialBarrier As List(Of RectangleF)
        Property DoubleSlit As RectangleF
        Property QuantumTunnelingPoints As List(Of PointF)
        Property InterferencePattern As List(Of Double)

        Sub New()
            PotentialBarrier = New List(Of RectangleF)
            QuantumTunnelingPoints = New List(Of PointF)
            InterferencePattern = New List(Of Double)
            
            ' Create a double-slit
            DoubleSlit = New RectangleF(400, 200, 20, 100)
            
            ' Create potential barriers
            PotentialBarrier.Add(New RectangleF(300, 150, 200, 10))
            PotentialBarrier.Add(New RectangleF(300, 350, 200, 10))
        End Sub

        Function CalculatePotentialAt(x As Single, y As Single) As Double
            ' Simple harmonic oscillator potential + barriers
            Dim centerX As Single = 400
            Dim centerY As Single = 300
            Dim distance As Double = Math.Sqrt((x - centerX) ^ 2 + (y - centerY) ^ 2)
            
            ' Harmonic potential: V = 0.5 * k * r^2
            Dim harmonic As Double = 0.0001 * distance * distance
            
            ' Add barrier potentials
            Dim barrierPotential As Double = 0
            For Each barrier In PotentialBarrier
                If x >= barrier.X AndAlso x <= barrier.X + barrier.Width AndAlso
                   y >= barrier.Y AndAlso y <= barrier.Y + barrier.Height Then
                    barrierPotential = 0.5
                    Exit For
                End If
            Next
            
            Return harmonic + barrierPotential
        End Function

        Function QuantumTunnel(particle As QuantumParticle, barrierHeight As Double) As Boolean
            ' Simple tunneling probability: P ≈ e^(-2 * κ * d)
            ' where κ = √(2m(V-E))/ħ
            Dim mass As Double = particle.Mass
            Dim energy As Double = 0.5 * mass * particle.Momentum.LengthSquared()
            
            If energy >= barrierHeight Then Return True ' Classical transmission
            
            Dim kappa As Double = Math.Sqrt(2 * mass * (barrierHeight - energy)) / PLANCK_CONSTANT
            Dim barrierWidth As Double = 10 ' Approximate width
            Dim tunnelingProb As Double = Math.Exp(-2 * kappa * barrierWidth)
            
            Return rnd.NextDouble() < tunnelingProb
        End Function
    End Class

    Class QuantumObserver
        Property Position As PointF
        Property MeasurementBasis As String ' "position", "momentum", "spin"
        Property MeasurementPrecision As Double ' 0 to 1
        Property IsMeasuring As Boolean
        Property DecoherenceRate As Double
        Property MeasurementAngle As Double ' For Stern-Gerlach

        Sub New()
            Position = New PointF(400, 300)
            MeasurementBasis = "position"
            MeasurementPrecision = 0.8
            DecoherenceRate = 0.1
            MeasurementAngle = 0
        End Sub

        Sub MeasureParticle(particle As QuantumParticle)
            ' Wave function collapse
            Dim measurementOutcome As Double
            
            Select Case MeasurementBasis
                Case "position"
                    ' Position measurement with uncertainty
                    Dim uncertainty As Double = (1 - MeasurementPrecision) * 50
                    measurementOutcome = particle.Position.X + (rnd.NextDouble() - 0.5) * uncertainty
                    particle.Position = New Vector2(CSng(measurementOutcome), particle.Position.Y)
                    
                Case "momentum"
                    ' Momentum measurement
                    measurementOutcome = particle.Momentum.Length()
                    particle.Momentum = Vector2.Normalize(particle.Momentum) * 
                                        CSng(measurementOutcome * MeasurementPrecision)
                    
                Case "spin"
                    ' Stern-Gerlach measurement
                    Dim angle As Double = MeasurementAngle * Math.PI / 180
                    Dim expectedSpin As Double = particle.Spin * Math.Cos(angle)
                    
                    ' Quantum probability: P(+) = cos²(θ/2), P(-) = sin²(θ/2)
                    Dim probUp As Double = Math.Pow(Math.Cos(angle / 2), 2)
                    particle.Spin = If(rnd.NextDouble() < probUp, 0.5, -0.5)
                    measurementOutcome = particle.Spin
            End Select
            
            ' Record measurement
            measurementHistory.Add(New MeasurementEvent(particle, MeasurementBasis, 
                                                       measurementOutcome, simulationTime))
            
            ' Decoherence
            particle.DecoherenceTime -= DecoherenceRate
        End Sub

        Function CalculateExpectationValue(particles As List(Of QuantumParticle), 
                                                 observable As String) As Double
            ' Calculate quantum expectation value
            Dim sum As Double = 0
            Dim count As Integer = 0
            
            For Each p In particles
                Select Case observable
                    Case "position"
                        sum += p.Position.X
                    Case "momentum"
                        sum += p.Momentum.Length()
                    Case "spin"
                        sum += p.Spin
                End Select
                count += 1
            Next
            
            Return If(count > 0, sum / count, 0)
        End Function
    End Class

    Class MeasurementEvent
        Property Particle As QuantumParticle
        Property Basis As String
        Property Value As Double
        Property Time As Double
        Property CollapsedState As Integer

        Sub New(p As QuantumParticle, b As String, v As Double, t As Double)
            Particle = p
            Basis = b
            Value = v
            Time = t
            CollapsedState = p.QuantumState
        End Sub
    End Class

    Declare Sub QuantumPhysicsForm_Load(sender As Object, e As EventArgs) Handles MyBase.Load
        DoubleBuffered = True
        InitializeQuantumSystem()
        
        simulationTimer.Interval = 16
        simulationTimer.Start()
    End Sub

    Declare Sub InitializeQuantumSystem()
        quantumParticles.Clear()
        entangledPairs.Clear()
        superpositionParticles.Clear()
        measurementHistory.Clear()

        ' Create quantum particles with different properties
        For i As Integer = 0 To 20
            Dim pos As New Vector2(rnd.Next(100, 700), rnd.Next(100, 500))
            Dim mom As New Vector2(CSng((rnd.NextDouble() - 0.5) * 2),
                                  CSng((rnd.NextDouble() - 0.5) * 2))
            Dim spin As Double = If(rnd.NextDouble() > 0.5, 0.5, -0.5)
            Dim charge As Double = If(rnd.NextDouble() > 0.5, 1, -1)
            
            Dim qp As New QuantumParticle(pos, mom, spin, charge)
            qp.Color = Color.FromArgb(200,
                CInt(128 + spin * 127),
                CInt(128 + charge * 127),
                255)
            
            quantumParticles.Add(qp)
            
            ' Create superposition for some particles
            If i < 5 Then
                superpositionParticles.Add(New SuperpositionState(qp))
            End If
        Next

        ' Create entangled pairs
        For i As Integer = 0 To 3
            Dim p1 As QuantumParticle = quantumParticles(i * 2)
            Dim p2 As QuantumParticle = quantumParticles(i * 2 + 1)
            entangledPairs.Add(New EntangledPair(p1, p2, "spin"))
        Next

        ' Set observer position
        observer.Position = New PointF(400, 300)
    End Sub

    Declare Sub UpdateQuantumSystem()
        simulationTime += 0.016

        ' Update particle wave functions
        For Each particle In quantumParticles
            ' Apply Schrödinger equation (simplified)
            Dim potential As Double = waveFunction.CalculatePotentialAt(particle.Position.X, particle.Position.Y)
            
            ' Time evolution: iħ ∂ψ/∂t = Hψ
            Dim kineticEnergy As Double = 0.5 * particle.Mass * particle.Momentum.LengthSquared()
            Dim totalEnergy As Complex = kineticEnergy + potential
            
            ' ψ(t+Δt) = ψ(t) * e^(-iEΔt/ħ)
            Dim phase As Double = -totalEnergy.Real * 0.016 / PLANCK_CONSTANT
            particle.WaveFunction *= New Complex(Math.Cos(phase), Math.Sin(phase))
            
            ' Update position based on probability current
            Dim probCurrent As Vector2 = particle.Momentum / CSng(particle.Mass)
            particle.Position += probCurrent * 0.016F
            
            ' Apply uncertainty principle occasionally
            If rnd.NextDouble() < 0.01 Then
                particle.ApplyUncertaintyPrinciple(5.0)
            End If
            
            ' Quantum tunneling check
            For Each barrier In waveFunction.PotentialBarrier
                If particle.Position.X >= barrier.X AndAlso particle.Position.X <= barrier.X + barrier.Width AndAlso
                   particle.Position.Y >= barrier.Y AndAlso particle.Position.Y <= barrier.Y + barrier.Height Then
                    
                    If waveFunction.QuantumTunnel(particle, 0.5) Then
                        ' Tunnel through
                        particle.Position = New Vector2(
                            particle.Position.X + barrier.Width * If(particle.Momentum.X > 0, 1, -1),
                            particle.Position.Y)
                        waveFunction.QuantumTunnelingPoints.Add(New PointF(particle.Position.X, particle.Position.Y))
                    Else
                        ' Reflect
                        particle.Momentum = New Vector2(-particle.Momentum.X, particle.Momentum.Y)
                    End If
                End If
            Next
            
            ' Double-slit interference
            If particle.Position.X >= waveFunction.DoubleSlit.X AndAlso
               particle.Position.X <= waveFunction.DoubleSlit.X + waveFunction.DoubleSlit.Width AndAlso
               particle.Position.Y >= waveFunction.DoubleSlit.Y AndAlso
               particle.Position.Y <= waveFunction.DoubleSlit.Y + waveFunction.DoubleSlit.Height Then
                
                ' Create interference pattern
                Dim slitWidth As Single = waveFunction.DoubleSlit.Width / 3
                Dim centerY As Single = waveFunction.DoubleSlit.Y + waveFunction.DoubleSlit.Height / 2
                
                If Math.Abs(particle.Position.Y - centerY) < slitWidth Then
                    ' Pass through central slit
                    ' Continue unchanged
                ElseIf particle.Position.Y < centerY - slitWidth AndAlso
                       particle.Position.Y > centerY - slitWidth * 3 Then
                    ' Pass through upper slit - add vertical momentum
                    particle.Momentum.Y += 0.5F
                ElseIf particle.Position.Y > centerY + slitWidth AndAlso
                       particle.Position.Y < centerY + slitWidth * 3 Then
                    ' Pass through lower slit - subtract vertical momentum
                    particle.Momentum.Y -= 0.5F
                Else
                    ' Hit the barrier - reflect
                    particle.Momentum = New Vector2(-particle.Momentum.X, particle.Momentum.Y)
                End If
            End If
        Next

        ' Update entanglement correlations
        For Each pair In entangledPairs
            pair.UpdateCorrelation()
        Next

        ' Decoherence
        For Each particle In quantumParticles
            particle.DecoherenceTime -= 0.001
            If particle.DecoherenceTime <= 0 Then
                particle.IsEntangled = False
                particle.QuantumState = If(rnd.NextDouble() > 0.5, 1, 0)
            End If
        Next

        ' Observer measurement if active
        If observer.IsMeasuring Then
            For Each particle In quantumParticles
                If Vector2.Distance(New Vector2(observer.Position.X, observer.Position.Y), 
                                   particle.Position) < 50 Then
                    observer.MeasureParticle(particle)
                End If
            Next
        End If

        ' Collapse superpositions when measured
        If isCollapsing Then
            For Each sup In superpositionParticles
                If Not sup.IsCollapsed Then
                    ' Collapse to one state based on probability
                    Dim rand As Double = rnd.NextDouble()
                    Dim cumulative As Double = 0
                    
                    For i As Integer = 0 To sup.ProbabilityDistribution.Count - 1
                        cumulative += sup.ProbabilityDistribution(i)
                        If rand <= cumulative Then
                            sup.BaseParticle.Position = sup.PossibleStates(i)
                            sup.IsCollapsed = True
                            
                            measurementHistory.Add(New MeasurementEvent(
                                sup.BaseParticle, "superposition", i, simulationTime))
                            Exit For
                        End If
                    Next
                End If
            Next
            isCollapsing = False
        End If

        ' Update interference pattern
        UpdateInterferencePattern()
    End Sub

    Declare Sub UpdateInterferencePattern()
        ' Simulate interference pattern on detector screen
        Dim detectorX As Integer = 700
        Dim intensities(200) As Double
        
        For Each particle In quantumParticles
            If particle.Position.X > 600 Then ' Particles past double slit
                Dim screenY As Integer = CInt(particle.Position.Y)
                If screenY >= 0 AndAlso screenY < intensities.Length Then
                    ' Calculate interference from two slits
                    Dim slit1Dist As Double = Math.Abs(screenY - 250) ' Upper slit
                    Dim slit2Dist As Double = Math.Abs(screenY - 350) ' Lower slit
                    
                    ' Path difference
                    Dim pathDiff As Double = Math.Abs(slit1Dist - slit2Dist)
                    Dim wavelength As Double = particle.CalculateDeBroglieWavelength()
                    
                    If wavelength > 0 Then
                        Dim phaseDiff As Double = 2 * Math.PI * pathDiff / wavelength
                        intensities(screenY) += Math.Pow(Math.Cos(phaseDiff / 2), 2)
                    End If
                End If
            End If
        Next
        
        waveFunction.InterferencePattern = intensities.ToList()
    End Sub

    Declare Sub QuantumPhysicsForm_Paint(sender As Object, e As PaintEventArgs) Handles MyBase.Paint
        e.Graphics.SmoothingMode = SmoothingMode.AntiAlias
        e.Graphics.Clear(Color.Black)

        ' Draw quantum information
        DrawQuantumInformation(e.Graphics)

        ' Draw potential barriers
        Using barrierPen As New Pen(Color.FromArgb(100, 255, 255, 0), 3)
            For Each barrier In waveFunction.PotentialBarrier
                e.Graphics.DrawRectangle(barrierPen, barrier.X, barrier.Y, barrier.Width, barrier.Height)
                e.Graphics.FillRectangle(New SolidBrush(Color.FromArgb(30, 255, 255, 0)),
                                        barrier.X, barrier.Y, barrier.Width, barrier.Height)
            Next
        End Using

        ' Draw double slit
        Using slitBrush As New SolidBrush(Color.FromArgb(100, 200, 200, 200))
            e.Graphics.FillRectangle(slitBrush, waveFunction.DoubleSlit)
            
            ' Draw slits
            Dim slitWidth As Single = waveFunction.DoubleSlit.Width / 3
            Dim centerY As Single = waveFunction.DoubleSlit.Y + waveFunction.DoubleSlit.Height / 2
            
            ' Upper slit
            e.Graphics.FillRectangle(Brushes.Black,
                                   waveFunction.DoubleSlit.X,
                                   centerY - slitWidth * 3,
                                   waveFunction.DoubleSlit.Width,
                                   slitWidth * 2)
            
            ' Lower slit
            e.Graphics.FillRectangle(Brushes.Black,
                                   waveFunction.DoubleSlit.X,
                                   centerY + slitWidth,
                                   waveFunction.DoubleSlit.Width,
                                   slitWidth * 2)
        End Using

        ' Draw interference pattern
        DrawInterferencePattern(e.Graphics)

        ' Draw superposition clouds
        For Each sup In superpositionParticles
            If Not sup.IsCollapsed Then
                For i As Integer = 0 To sup.PossibleStates.Count - 1
                    Dim alpha As Integer = CInt(sup.ProbabilityDistribution(i) * 100)
                    Using brush As New SolidBrush(Color.FromArgb(alpha, 100, 200, 255))
                        e.Graphics.FillEllipse(brush,
                                              sup.PossibleStates(i).X - 10,
                                              sup.PossibleStates(i).Y - 10,
                                              20, 20)
                    End Using
                Next
            End If
        Next

        ' Draw particles
        For Each particle In quantumParticles
            ' Draw wave function magnitude as glow
            Dim waveMagnitude As Double = particle.WaveFunction.Magnitude
            Dim glowSize As Single = particle.Size + CSng(waveMagnitude * 10)
            
            Using glowBrush As New SolidBrush(Color.FromArgb(50, particle.Color))
                e.Graphics.FillEllipse(glowBrush,
                                      particle.Position.X - glowSize / 2,
                                      particle.Position.Y - glowSize / 2,
                                      glowSize, glowSize)
            End Using

            ' Draw particle
            Using particleBrush As New SolidBrush(particle.Color)
                e.Graphics.FillEllipse(particleBrush,
                                      particle.Position.X - particle.Size / 2,
                                      particle.Position.Y - particle.Size / 2,
                                      particle.Size, particle.Size)
            End Using

            ' Draw spin indicator
            If Math.Abs(particle.Spin) > 0 Then
                Dim spinColor As Color = If(particle.Spin > 0, Color.Cyan, Color.Magenta)
                Using spinPen As New Pen(spinColor, 2)
                    e.Graphics.DrawLine(spinPen,
                                       particle.Position.X,
                                       particle.Position.Y - 5,
                                       particle.Position.X,
                                       particle.Position.Y + 5)
                End Using
            End If
        Next

        ' Draw entanglement connections
        For Each pair In entangledPairs
            Using entanglePen As New Pen(Color.FromArgb(100, 0, 255, 0), 2)
                e.Graphics.DrawLine(entanglePen,
                                   pair.ParticleA.Position.X,
                                   pair.ParticleA.Position.Y,
                                   pair.ParticleB.Position.X,
                                   pair.ParticleB.Position.Y)
            End Using
        End Next

        ' Draw observer
        DrawObserver(e.Graphics)

        ' Draw measurement history
        DrawMeasurements(e.Graphics)

        ' Draw quantum tunneling points
        For Each tunnelPoint In waveFunction.QuantumTunnelingPoints
            Using tunnelBrush As New SolidBrush(Color.FromArgb(100, 255, 165, 0))
                e.Graphics.FillEllipse(tunnelBrush,
                                      tunnelPoint.X - 3,
                                      tunnelPoint.Y - 3,
                                      6, 6)
            End Using
        Next
    End Sub

    Declare Sub DrawQuantumInformation(g As Graphics)
        ' Draw quantum theory information
        Using font As New Font("Consolas", 9)
            Dim yPos As Integer = 10
            
            g.DrawString($"QuantumLoop Simulation", New Font("Arial", 12, FontStyle.Bold), 
                        Brushes.White, 10, yPos)
            yPos += 25
            
            g.DrawString($"Time: {simulationTime:F2}s | Particles: {quantumParticles.Count}", 
                        font, Brushes.White, 10, yPos)
            yPos += 20
            
            g.DrawString($"Superpositions: {superpositionParticles.Count} | Entangled Pairs: {entangledPairs.Count}", 
                        font, Brushes.Cyan, 10, yPos)
            yPos += 20
            
            g.DrawString($"Measurement Basis: {observer.MeasurementBasis} | Precision: {observer.MeasurementPrecision:P0}", 
                        font, Brushes.Yellow, 10, yPos)
            yPos += 20
            
            ' Calculate quantum properties
            Dim avgSpin As Double = observer.CalculateExpectationValue(quantumParticles, "spin")
            Dim avgMomentum As Double = observer.CalculateExpectationValue(quantumParticles, "momentum")
            
            g.DrawString($"⟨S⟩ = {avgSpin:F3} | ⟨p⟩ = {avgMomentum:E2} kg⋅m/s", 
                        font, Brushes.LightGreen, 10, yPos)
            yPos += 20
            
            ' Quantum principles
            g.DrawString("Quantum Principles:", New Font("Arial", 10, FontStyle.Bold), 
                        Brushes.White, 10, yPos)
            yPos += 20
            
            Dim principles As String() = {
                "• Superposition: Particles exist in multiple states simultaneously",
                "• Entanglement: Instant correlation between distant particles",
                "• Wave Function Collapse: Measurement determines reality",
                "• Quantum Tunneling: Particles penetrate classically forbidden barriers",
                "• Uncertainty Principle: Δx⋅Δp ≥ ħ/2",
                "• Decoherence: Quantum to classical transition"
            }
            
            For Each principle In principles
                g.DrawString(principle, font, Brushes.LightGray, 20, yPos)
                yPos += 16
            Next
        End Using
    End Sub

    Declare Sub DrawObserver(g As Graphics)
        ' Draw observer with measurement cone
        Using observerBrush As New SolidBrush(Color.FromArgb(200, 255, 255, 255))
            g.FillEllipse(observerBrush,
                         observer.Position.X - 10,
                         observer.Position.Y - 10,
                         20, 20)
        End Using
        
        g.DrawString("OBSERVER", New Font("Arial", 8), Brushes.White,
                    observer.Position.X - 25, observer.Position.Y + 15)

        ' Draw measurement cone if measuring
        If observer.IsMeasuring Then
            Using measurementBrush As New SolidBrush(Color.FromArgb(50, 255, 255, 0))
                g.FillPie(measurementBrush,
                         observer.Position.X - 50,
                         observer.Position.Y - 50,
                         100, 100,
                         CSng(observer.MeasurementAngle - 30),
                         60)
            End Using
        End If
    End Sub

    Declare Sub DrawInterferencePattern(g As Graphics)
        ' Draw detector screen
        Dim screenX As Integer = 700
        g.DrawLine(Pens.Gray, screenX, 0, screenX, 600)
        
        ' Draw interference pattern
        If waveFunction.InterferencePattern.Count > 0 Then
            For y As Integer = 0 To waveFunction.InterferencePattern.Count - 1
                Dim intensity As Double = waveFunction.InterferencePattern(y)
                If intensity > 0 Then
                    Dim brightness As Integer = CInt(Math.Min(intensity * 50, 255))
                    Using patternPen As New Pen(Color.FromArgb(brightness, 255, 255, 255), 2)
                        g.DrawLine(patternPen,
                                  screenX,
                                  y,
                                  screenX + CInt(intensity * 20),
                                  y)
                    End Using
                End If
            Next
        End If
        
        g.DrawString("INTERFERENCE PATTERN", New Font("Arial", 8), Brushes.White, screenX + 5, 10)
    End Sub

    Declare Sub DrawMeasurements(g As Graphics)
        ' Draw recent measurement events
        Dim recentMeasurements = measurementHistory.
            OrderByDescending(Function(m) m.Time).
            Take(5).ToList()
        
        Using font As New Font("Consolas", 8)
            For i As Integer = 0 To recentMeasurements.Count - 1
                Dim m As MeasurementEvent = recentMeasurements(i)
                Dim text As String = $"{m.Basis}: {m.Value:F3} at t={m.Time:F2}s"
                g.DrawString(text, font, Brushes.Yellow, 600, 400 + i * 15)
            Next
        End Using
    End Sub

    Declare Sub simulationTimer_Tick(sender As Object, e As EventArgs) Handles simulationTimer.Tick
        UpdateQuantumSystem()
        Invalidate()
    End Sub

    ' Control buttons
    Declare Sub btnCollapse_Click(sender As Object, e As EventArgs) Handles btnCollapse.Click
        isCollapsing = True
        For Each sup In superpositionParticles
            sup.IsCollapsed = False
        Next
    End Sub

    Declare Sub btnEntangle_Click(sender As Object, e As EventArgs) Handles btnEntangle.Click
        ' Create new entangled pairs
        If quantumParticles.Count >= 2 Then
            For i As Integer = 0 To Math.Min(3, quantumParticles.Count \ 2) - 1
                Dim p1 As QuantumParticle = quantumParticles(i * 2)
                Dim p2 As QuantumParticle = quantumParticles(i * 2 + 1)
                entangledPairs.Add(New EntangledPair(p1, p2, "spin"))
            Next
        End If
    End Sub

    Declare Sub btnMeasure_Click(sender As Object, e As EventArgs) Handles btnMeasure.Click
        observer.IsMeasuring = Not observer.IsMeasuring
        btnMeasure.Text = If(observer.IsMeasuring, "Stop Measuring", "Start Measuring")
    End Sub

    Declare Sub btnChangeBasis_Click(sender As Object, e As EventArgs) Handles btnChangeBasis.Click
        Dim bases As String() = {"position", "momentum", "spin"}
        Dim currentIndex As Integer = Array.IndexOf(bases, observer.MeasurementBasis)
        observer.MeasurementBasis = bases((currentIndex + 1) Mod bases.Length)
        
        If observer.MeasurementBasis = "spin" Then
            observer.MeasurementAngle = rnd.Next(0, 360)
        End If
    End Sub

    Declare Sub btnReset_Click(sender As Object, e As EventArgs) Handles btnReset.Click
        InitializeQuantumSystem()
        simulationTime = 0
    End Sub

    Declare Sub QuantumPhysicsForm_KeyDown(sender As Object, e As KeyEventArgs) Handles MyBase.KeyDown
        Select Case e.KeyCode
            Case Keys.Space
                ' Toggle measurement
                btnMeasure_Click(Nothing, Nothing)
            Case Keys.C
                ' Collapse wave function
                btnCollapse_Click(Nothing, Nothing)
            Case Keys.E
                ' Create entanglement
                btnEntangle_Click(Nothing, Nothing)
            Case Keys.B
                ' Change measurement basis
                btnChangeBasis_Click(Nothing, Nothing)
            Case Keys.R
                ' Reset simulation
                btnReset_Click(Nothing, Nothing)
            Case Keys.Add, Keys.Oemplus
                ' Increase measurement precision
                observer.MeasurementPrecision = Math.Min(observer.MeasurementPrecision + 0.1, 1.0)
            Case Keys.Subtract, Keys.OemMinus
                ' Decrease measurement precision
                observer.MeasurementPrecision = Math.Max(observer.MeasurementPrecision - 0.1, 0.1)
        End Select
    End Sub

End Class
