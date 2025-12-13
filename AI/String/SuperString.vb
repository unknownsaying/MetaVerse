Imports System.Drawing
Imports System.Drawing.Drawing2D

Public Class StringCycle
    Private stringsList As New List(Of GuitarString)
    Private timer As New Timer()
    Private time As Double = 0
    Private amplitudeFactor As Double = 30
    Private stringTension As Double = 1.0
    Private supersymmetryMode As Boolean = False
    Private extraDimensions As Integer = 6

    Public Sub New()
        InitializeComponent()
        InitializeStrings()
        SetupTimer()
        SetupUI()
    End Sub

    Private Sub InitializeStrings()
        ' Standard tuning frequencies (E2, A2, D3, G3, B3, E4)
        Dim frequencies() As Double = {82.41, 110.0, 146.83, 196.0, 246.94, 329.63}
        Dim notes() As String = {"E", "A", "D", "G", "B", "E"}
        Dim colors() As Color = {Color.Red, Color.Orange, Color.Yellow, 
                                 Color.Green, Color.Blue, Color.Purple}
        
        For i As Integer = 0 To 5
            stringsList.Add(New GuitarString With {
                .Note = notes(i),
                .Frequency = frequencies(i),
                .StringColor = colors(i),
                .YPosition = 100 + i * 60,
                .Length = 700,
                .Tension = stringTension,
                .MassDensity = 0.001 * (6 - i) ' Higher strings have less mass density
            })
        Next
    End Sub

    Private Sub SetupTimer()
        timer.Interval = 16 ' ~60 FPS
        AddHandler timer.Tick, AddressOf Timer_Tick
        timer.Start()
    End Sub

    Private Sub SetupUI()
        Me.Text = "Superstring Guitar Theory Visualizer"
        Me.Size = New Size(900, 700)
        Me.DoubleBuffered = True
        
        ' Add controls
        AddHandler Me.Paint, AddressOf MainForm_Paint
        AddHandler Me.KeyDown, AddressOf MainForm_KeyDown
        
        ' Add buttons
        Dim pluckBtn As New Button With {
            .Text = "Pluck All Strings",
            .Location = New Point(10, 10),
            .Size = New Size(120, 30)
        }
        AddHandler pluckBtn.Click, AddressOf PluckStrings
        Me.Controls.Add(pluckBtn)
        
        Dim supersymmetryBtn As New Button With {
            .Text = "Toggle Supersymmetry",
            .Location = New Point(140, 10),
            .Size = New Size(150, 30)
        }
        AddHandler supersymmetryBtn.Click, AddressOf ToggleSupersymmetry
        Me.Controls.Add(supersymmetryBtn)
        
        Dim tensionTrack As New TrackBar With {
            .Location = New Point(300, 10),
            .Size = New Size(200, 30),
            .Minimum = 1,
            .Maximum = 100,
            .Value = CInt(stringTension * 10)
        }
        AddHandler tensionTrack.ValueChanged, AddressOf TensionChanged
        Me.Controls.Add(tensionTrack)
        
        AddHandler Me.Resize, AddressOf MainForm_Resize
    End Sub

    Private Sub Timer_Tick(sender As Object, e As EventArgs)
        time += 0.05
        Me.Invalidate()
    End Sub

    Private Sub MainForm_Paint(sender As Object, e As PaintEventArgs)
        e.Graphics.SmoothingMode = SmoothingMode.AntiAlias
        e.Graphics.Clear(Color.Black)
        
        ' Draw superstring dimensions grid
        DrawDimensionGrid(e.Graphics)
        
        ' Draw strings with vibrations
        For Each gs In stringsList
            DrawString(e.Graphics, gs)
        Next
        
        ' Draw string information
        DrawInfo(e.Graphics)
        
        ' If supersymmetry mode is on, draw partner strings
        If supersymmetryMode Then
            DrawSupersymmetryStrings(e.Graphics)
        End If
    End Sub

    Private Sub DrawDimensionGrid(g As Graphics)
        Dim pen As New Pen(Color.FromArgb(30, 100, 100, 100), 1)
        
        ' Draw Calabi-Yau manifold inspired patterns
        For i As Integer = 0 To 360 Step 30
            Dim angle As Double = i * Math.PI / 180
            Dim radius As Integer = 50
            Dim center As New Point(Me.ClientSize.Width - 150, 150)
            
            For j As Integer = 1 To 3
                Dim r As Integer = radius * j
                Dim x1 As Integer = center.X + CInt(r * Math.Cos(angle))
                Dim y1 As Integer = center.Y + CInt(r * Math.Sin(angle))
                
                ' Draw dimensional circles
                g.DrawEllipse(pen, center.X - r, center.Y - r, 2 * r, 2 * r)
            Next
        Next
    End Sub

    Private Sub DrawString(g As Graphics, gs As GuitarString)
        Dim points As New List(Of Point)
        Dim segments As Integer = 100
        
        ' Calculate wave equation based on string physics
        For i As Integer = 0 To segments
            Dim x As Double = 100 + (gs.Length * i / segments)
            Dim t As Double = time
            
            ' Wave equation for vibrating string
            Dim yOffset As Double = CalculateStringWave(i / CDbl(segments), t, gs)
            
            ' Add extra dimension effects
            If extraDimensions > 0 Then
                yOffset += 5 * Math.Sin(2 * Math.PI * gs.Frequency * t / 100 + 
                                      i * Math.PI / segments * extraDimensions)
            End If
            
            points.Add(New Point(CInt(x), CInt(gs.YPosition + yOffset)))
        Next
        
        ' Draw the string
        Dim stringPen As New Pen(gs.StringColor, 2)
        g.DrawLines(stringPen, points.ToArray())
        
        ' Draw vibration particles (string theory particles)
        DrawVibrationParticles(g, points)
    End Sub

    Private Function CalculateStringWave(x As Double, t As Double, gs As GuitarString) As Double
        ' Standing wave equation for a vibrating string
        Dim k As Double = Math.PI / (gs.Length / 100) ' Wave number
        Dim omega As Double = 2 * Math.PI * gs.Frequency ' Angular frequency
        
        ' Damped oscillation to simulate energy loss
        Dim damping As Double = Math.Exp(-0.01 * t)
        
        ' Fundamental mode + harmonics
        Dim y As Double = 0
        For harmonic As Integer = 1 To 5
            y += (amplitudeFactor / harmonic) * 
                 Math.Sin(harmonic * k * x * 10) * 
                 Math.Cos(harmonic * omega * t / 100) *
                 damping
        Next
        
        Return y * (gs.Tension / gs.MassDensity) * 0.1
    End Function

    Private Sub DrawVibrationParticles(g As Graphics, points As List(Of Point))
        Dim rnd As New Random()
        
        ' Draw particles at antinodes
        For i As Integer = 10 To points.Count - 1 Step 20
            If rnd.NextDouble() > 0.7 Then
                Dim brush As New SolidBrush(Color.FromArgb(150, Color.White))
                Dim size As Integer = 3 + rnd.Next(5)
                g.FillEllipse(brush, points(i).X - size \ 2, 
                                       points(i).Y - size \ 2, size, size)
            End If
        Next
    End Sub

    Private Sub DrawSupersymmetryStrings(g As Graphics)
        ' Draw supersymmetric partner strings
        Dim pen As New Pen(Color.FromArgb(100, Color.Magenta), 1)
        pen.DashStyle = DashStyle.Dash
        
        For Each gs In stringsList
            Dim points As New List(Of Point)
            Dim segments As Integer = 100
            
            For i As Integer = 0 To segments
                Dim x As Double = 100 + (gs.Length * i / segments)
                Dim t As Double = time
                
                ' Supersymmetric partner vibrates out of phase
                Dim yOffset As Double = -CalculateStringWave(i / CDbl(segments), t, gs) * 0.5
                yOffset += 10 * Math.Sin(2 * Math.PI * gs.Frequency * t / 50)
                
                points.Add(New Point(CInt(x), CInt(gs.YPosition + yOffset)))
            Next
            
            g.DrawLines(pen, points.ToArray())
        Next
    End Sub
    Private Sub DrawInfo(g As Graphics)
        Dim font As New Font("Consolas", 10)
        Dim brush As New SolidBrush(Color.White)
        
        ' String theory info
        g.DrawString("Superstring Theory Visualization", font, brush, 10, 50)
        g.DrawString($"Extra Dimensions: {extraDimensions}", font, brush, 10, 70)
        g.DrawString($"Tension: {stringTension:F2}", font, brush, 10, 90)
        
        ' Guitar tuning info
        For i As Integer = 0 To stringsList.Count - 1
            Dim gs = stringsList(i)
            Dim info As String = $"String {6 - i}: {gs.Note} ({gs.Frequency:F2} Hz)"
            g.DrawString(info, font, New SolidBrush(gs.StringColor), 10, 120 + i * 20)
        Next
        
        ' Physics info
        g.DrawString("Wave Equation: ∂²y/∂t² = (T/μ) ∂²y/∂x²", font, brush, 500, 50)
        g.DrawString("T = Tension, μ = Mass per unit length", font, brush, 500, 70)
    End Sub

    ' Event Handlers
    Private Sub PluckStrings(sender As Object, e As EventArgs)
        amplitudeFactor = 50
        time = 0
        For Each gs In stringsList
            gs.Tension = stringTension * (1 + 0.1 * (New Random()).NextDouble())
        Next
    End Sub

    Private Sub ToggleSupersymmetry(sender As Object, e As EventArgs)
        supersymmetryMode = Not supersymmetryMode
        extraDimensions = If(supersymmetryMode, 10, 6)
    End Sub

    Private Sub TensionChanged(sender As Object, e As EventArgs)
        Dim trackBar = DirectCast(sender, TrackBar)
        stringTension = trackBar.Value / 10.0
        For Each gs In stringsList
            gs.Tension = stringTension
        Next
    End Sub

    Private Sub MainForm_KeyDown(sender As Object, e As KeyEventArgs)
        Select Case e.KeyCode
            Case Keys.Space
                PluckStrings(Nothing, Nothing)
            Case Keys.S
                ToggleSupersymmetry(Nothing, Nothing)
            Case Keys.Add
                extraDimensions = Math.Min(11, extraDimensions + 1)
            Case Keys.Subtract
                extraDimensions = Math.Max(0, extraDimensions - 1)
        End Select
    End Sub

    Private Sub MainForm_Resize(sender As Object, e As EventArgs)
        Me.Invalidate()
    End Sub
End Class