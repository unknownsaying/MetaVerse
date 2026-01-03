Imports System.Drawing.Drawing2D
Imports System.Math

Class MainForm
    simulationType As String = "Loop"
    simulationRunning As Boolean = False
    simulationTime As Double = 0
    particles As List(Of Particle) = New List(Of Particle)()
    rnd As Random = New Random()

    ' Physics parameters
    gravity As Single = 0.1
    friction As Single = 0.99
    quantumProbability As Double = 0.05
    entanglementDistance As Integer = 100
    decayRate As Single = 0.001

    Class Particle
        Property Position As PointF
        Property Velocity As PointF
        Property Color As Color
        Property Size As Single
        Property Life As Single = 1.0
        Property QuantumState As Integer = 0
        Property EntangledWith As List(Of Particle) = New List(Of Particle)()
        Property IsObserver As Boolean = False
        Property DecayCounter As Single = 0

        Sub New(pos As PointF, vel As PointF, clr As Color, sz As Single)
            Position = pos
            Velocity = vel
            Color = clr
            Size = sz
        End Sub
    End Class

    Sub MainForm_Load(sender As Object, e As EventArgs) Handles MyBase.Load
        DoubleBuffered = True
        InitializeParticles()
        Timer1.Interval = 16 ' ~60 FPS
        Timer1.Start()
    End Sub

    Sub InitializeParticles()
        particles.Clear()
        Dim centerX As Integer = pbCanvas.Width \ 2
        Dim centerY As Integer = pbCanvas.Height \ 2

        Select Case simulationType
            Case "Loop"
                ' Create particles in a circular loop
                For i As Integer = 0 To 50
                    Dim angle As Single = CSng(i) / 50 * Math.PI * 2
                    Dim radius As Integer = 150
                    Dim pos As New PointF(
                        centerX + CSng(Math.Cos(angle)) * radius,
                        centerY + CSng(Math.Sin(angle)) * radius)
                    Dim vel As New PointF(
                        CSng(Math.Sin(angle)) * 0.5,
                        -CSng(Math.Cos(angle)) * 0.5)
                    particles.Add(New Particle(pos, vel, Color.FromArgb(150, 0, 150, 255), 4))
                Next

            Case "QuantumLoop"
                ' Create entangled particles
                For i As Integer = 0 To 30
                    Dim angle As Single = CSng(i) / 30 * Math.PI * 2
                    Dim radius As Integer = 100
                    Dim pos As New PointF(
                        centerX + CSng(Math.Cos(angle)) * radius,
                        centerY + CSng(Math.Sin(angle)) * radius)
                    Dim vel As New PointF(
                        CSng(Math.Sin(angle)) * 0.3,
                        -CSng(Math.Cos(angle)) * 0.3)
                    Dim p As New Particle(pos, vel, Color.FromArgb(150, 255, 0, 200), 5)
                    p.QuantumState = rnd.Next(0, 2)
                    particles.Add(p)
                Next

                ' Entangle particles
                For Each p1 In particles
                    For Each p2 In particles
                        If p1 IsNot p2 AndAlso Distance(p1.Position, p2.Position) < entanglementDistance Then
                            If Not p1.EntangledWith.Contains(p2) Then
                                p1.EntangledWith.Add(p2)
                                p2.EntangledWith.Add(p1)
                            End If
                        End If
                    Next
                Next

            Case "DeathLoop"
                ' Create decaying particles
                For i As Integer = 0 To 40
                    Dim angle As Single = CSng(i) / 40 * Math.PI * 2
                    Dim radius As Integer = 120
                    Dim pos As New PointF(
                        centerX + CSng(Math.Cos(angle)) * radius,
                        centerY + CSng(Math.Sin(angle)) * radius)
                    Dim vel As New PointF(
                        CSng(Math.Sin(angle)) * 0.4,
                        -CSng(Math.Cos(angle)) * 0.4)
                    Dim p As New Particle(pos, vel, Color.FromArgb(200, 255, 50, 50), 6)
                    p.DecayCounter = rnd.Next(100, 500)
                    particles.Add(p)
                Next
        End Select

        ' Add observer particle
        Dim observer As New Particle(
            New PointF(centerX, centerY - 200),
            New PointF(0, 0),
            Color.White, 8)
        observer.IsObserver = True
        particles.Add(observer)
    End Sub

    Function Distance(p1 As PointF, p2 As PointF) As Single
        Dim dx As Single = p1.X - p2.X
        Dim dy As Single = p1.Y - p2.Y
        Return CSng(Math.Sqrt(dx * dx + dy * dy))
    End Function

    Sub UpdatePhysics()
        simulationTime += 0.016

        Select Case simulationType
            Case "Loop"
                UpdateLoopPhysics()
            Case "QuantumLoop"
                UpdateQuantumPhysics()
            Case "DeathLoop"
                UpdateDeathLoopPhysics()
        End Select
    End Sub

    Sub UpdateLoopPhysics()
        For Each p In particles
            If Not p.IsObserver Then
                ' Apply centripetal force
                Dim center As New PointF(pbCanvas.Width \ 2, pbCanvas.Height \ 2)
                Dim dir As New PointF(p.Position.X - center.X, p.Position.Y - center.Y)
                Dim distance As Single = Math.Max(Distance(p.Position, center), 1)
                Dim forceMagnitude As Single = (distance - 150) * 0.01

                p.Velocity = New PointF(
                    (p.Velocity.X - dir.X / distance * forceMagnitude) * friction,
                    (p.Velocity.Y - dir.Y / distance * forceMagnitude) * friction)

                ' Add slight randomness
                p.Velocity = New PointF(
                    p.Velocity.X + (rnd.NextSingle() - 0.5) * 0.05,
                    p.Velocity.Y + (rnd.NextSingle() - 0.5) * 0.05)

                p.Position = New PointF(
                    p.Position.X + p.Velocity.X,
                    p.Position.Y + p.Velocity.Y)

                ' Color oscillation
                Dim hue As Integer = CInt((Math.Sin(simulationTime * 2 + p.Position.X * 0.01) + 1) * 127)
                p.Color = Color.FromArgb(150, hue, 150, 255)
            End If
        Next
    End Sub

    Sub UpdateQuantumPhysics()
        ' Quantum superposition and entanglement effects
        For Each p In particles
            If Not p.IsObserver Then
                ' Quantum state fluctuations
                If rnd.NextDouble() < quantumProbability Then
                    p.QuantumState = 1 - p.QuantumState ' Flip state
                    ' Entanglement effect
                    For Each entangled In p.EntangledWith
                        entangled.QuantumState = p.QuantumState
                    Next
                End If

                ' Wave function behavior
                Dim center As New PointF(pbCanvas.Width \ 2, pbCanvas.Height \ 2)
                Dim dir As New PointF(p.Position.X - center.X, p.Position.Y - center.Y)
                Dim distance As Single = Distance(p.Position, center)
                Dim targetRadius As Single = 100 + Math.Sin(simulationTime + p.QuantumState * Math.PI) * 30

                Dim forceMagnitude As Single = (distance - targetRadius) * 0.015

                p.Velocity = New PointF(
                    (p.Velocity.X - dir.X / Math.Max(distance, 1) * forceMagnitude) * friction,
                    (p.Velocity.Y - dir.Y / Math.Max(distance, 1) * forceMagnitude) * friction)

                p.Position = New PointF(
                    p.Position.X + p.Velocity.X,
                    p.Position.Y + p.Velocity.Y)

                ' Quantum color based on state
                p.Color = Color.FromArgb(150,
                    p.QuantumState * 255,
                    100,
                    (1 - p.QuantumState) * 200 + 55)
            End If
        Next
    End Sub

    Sub UpdateDeathLoopPhysics()
        For Each p In particles
            If Not p.IsObserver Then
                ' Decay process
                p.DecayCounter -= 1
                If p.DecayCounter <= 0 Then
                    p.Life -= decayRate
                    p.Size *= 0.995

                    ' Spawn new particles when decaying
                    If rnd.NextDouble() < 0.02 AndAlso p.Life > 0.1 Then
                        Dim newParticle As New Particle(
                            New PointF(p.Position.X, p.Position.Y),
                            New PointF((rnd.NextSingle() - 0.5) * 2, (rnd.NextSingle() - 0.5) * 2),
                            Color.FromArgb(100, 255, 100, 0),
                            3)
                        particles.Add(newParticle)
                    End If
                End If

                ' Chaotic orbit with decay
                Dim center As New PointF(pbCanvas.Width \ 2, pbCanvas.Height \ 2)
                Dim angle As Single = CSng(Math.Atan2(p.Position.Y - center.Y, p.Position.X - center.X))
                Dim radius As Single = Distance(p.Position, center)
                Dim targetRadius As Single = 120 + Math.Sin(simulationTime * 0.5 + p.Position.X * 0.01) * 50

                ' Apply decay-modified centripetal force
                Dim forceMagnitude As Single = (radius - targetRadius) * 0.01 * p.Life

                p.Velocity = New PointF(
                    (p.Velocity.X - (p.Position.X - center.X) / Math.Max(radius, 1) * forceMagnitude) * (friction * p.Life),
                    (p.Velocity.Y - (p.Position.Y - center.Y) / Math.Max(radius, 1) * forceMagnitude) * (friction * p.Life))

                ' Add chaos
                p.Velocity = New PointF(
                    p.Velocity.X + (rnd.NextSingle() - 0.5) * 0.1 * (1 - p.Life),
                    p.Velocity.Y + (rnd.NextSingle() - 0.5) * 0.1 * (1 - p.Life))

                p.Position = New PointF(
                    p.Position.X + p.Velocity.X,
                    p.Position.Y + p.Velocity.Y)

                ' Decay-based color
                p.Color = Color.FromArgb(CInt(p.Life * 200),
                    CInt(255 * (1 - p.Life) * 0.5 + 200 * p.Life),
                    CInt(50 * p.Life),
                    CInt(50 * p.Life))

                ' Remove dead particles
                If p.Life <= 0 Then
                    p.Life = 0
                    p.Size = 0
                End If
            End If
        Next

        ' Remove particles with zero size
        particles.RemoveAll(Function(x) x.Size <= 0)
    End Sub

    Sub pbCanvas_Paint(sender As Object, e As PaintEventArgs) Handles pbCanvas.Paint
        e.Graphics.SmoothingMode = SmoothingMode.AntiAlias
        e.Graphics.Clear(Color.Black)

        ' Draw time indicator
        Dim timeText As String = $"Time: {simulationTime:F2}s | Particles: {particles.Count}"
        e.Graphics.DrawString(timeText, Me.Font, Brushes.White, 10, 10)

        ' Draw physics theory labels
        Select Case simulationType
            Case "Loop"
                e.Graphics.DrawString("Classical Loop - Stable Orbital Mechanics", 
                                     New Font("Arial", 10), Brushes.Cyan, 10, 30)
                e.Graphics.DrawString("• Conservation of Angular Momentum", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 50)
                e.Graphics.DrawString("• Newtonian Gravity", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 65)
            Case "QuantumLoop"
                e.Graphics.DrawString("Quantum Loop - Superposition & Entanglement", 
                                     New Font("Arial", 10), Brushes.Magenta, 10, 30)
                e.Graphics.DrawString("• Quantum Superposition", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 50)
                e.Graphics.DrawString("• Quantum Entanglement", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 65)
                e.Graphics.DrawString("• Wave Function Collapse", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 80)
            Case "DeathLoop"
                e.Graphics.DrawString("Death Loop - Entropy & Decay", 
                                     New Font("Arial", 10), Brushes.OrangeRed, 10, 30)
                e.Graphics.DrawString("• Entropy Increase", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 50)
                e.Graphics.DrawString("• Particle Decay", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 65)
                e.Graphics.DrawString("• Chaotic Dynamics", 
                                     New Font("Arial", 8), Brushes.LightGray, 20, 80)
        End Select

        ' Draw particles
        For Each p In particles
            If p.Life > 0 Then
                Using brush As New SolidBrush(p.Color)
                    e.Graphics.FillEllipse(brush, 
                                          p.Position.X - p.Size / 2,
                                          p.Position.Y - p.Size / 2,
                                          p.Size, p.Size)

                    ' Draw observer differently
                    If p.IsObserver Then
                        e.Graphics.DrawEllipse(Pens.White, 
                                             p.Position.X - p.Size,
                                             p.Position.Y - p.Size,
                                             p.Size * 2, p.Size * 2)
                        e.Graphics.DrawString("Observer", 
                                             New Font("Arial", 8), 
                                             Brushes.White,
                                             p.Position.X - 20, p.Position.Y + p.Size + 5)
                    End If
                End Using
            End If
        Next

        ' Draw connections for QuantumLoop
        If simulationType = "QuantumLoop" Then
            For Each p In particles
                If Not p.IsObserver Then
                    For Each entangled In p.EntangledWith
                        If particles.IndexOf(entangled) > particles.IndexOf(p) Then
                            Using pen As New Pen(Color.FromArgb(50, 200, 100, 255), 1)
                                e.Graphics.DrawLine(pen, p.Position, entangled.Position)
                            End Using
                        End If
                    Next
                End If
            Next
        End If

        ' Draw center point
        Dim centerX As Integer = pbCanvas.Width \ 2
        Dim centerY As Integer = pbCanvas.Height \ 2
        e.Graphics.DrawEllipse(Pens.Gray, centerX - 3, centerY - 3, 6, 6)
    End Sub

    Sub Timer1_Tick(sender As Object, e As EventArgs) Handles Timer1.Tick
        If simulationRunning Then
            UpdatePhysics()
            pbCanvas.Invalidate()
        End If
    End Sub

    Sub btnLoop_Click(sender As Object, e As EventArgs) Handles btnLoop.Click
        simulationType = "Loop"
        simulationTime = 0
        InitializeParticles()
        pbCanvas.Invalidate()
    End Sub

    Sub btnQuantumLoop_Click(sender As Object, e As EventArgs) Handles btnQuantumLoop.Click
        simulationType = "QuantumLoop"
        simulationTime = 0
        InitializeParticles()
        pbCanvas.Invalidate()
    End Sub

    Sub btnDeathLoop_Click(sender As Object, e As EventArgs) Handles btnDeathLoop.Click
        simulationType = "DeathLoop"
        simulationTime = 0
        InitializeParticles()
        pbCanvas.Invalidate()
    End Sub

    Sub btnStartStop_Click(sender As Object, e As EventArgs) Handles btnStartStop.Click
        simulationRunning = Not simulationRunning
        btnStartStop.Text = If(simulationRunning, "Pause", "Start")
    End Sub

    Sub btnReset_Click(sender As Object, e As EventArgs) Handles btnReset.Click
        simulationTime = 0
        InitializeParticles()
        pbCanvas.Invalidate()
    End Sub

End Class
