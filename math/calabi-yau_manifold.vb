' CalabiYauVisualization.vb
' 6-Dimensional Calabi-Yau Manifold Compactification Visualizer
' String Theory & Supersymmetry Visualization

Imports System
Imports System.Drawing
Imports System.Drawing.Drawing2D
Imports System.Numerics
Imports System.Windows.Forms
Imports System.Collections.Generic
Imports System.Linq

' ==============================================
' 1. MATHEMATICAL FOUNDATIONS
' ==============================================

''' <summary>
''' Represents a 6-dimensional Calabi-Yau manifold
''' Compactification of extra dimensions in string theory
''' </summary>
Public Class CalabiYauManifold
    
    ' Calabi-Yau Properties
    Public Property Dimension As Integer = 6
    Public Property EulerCharacteristic As Integer
    Public Property HodgeDiamond As Integer(,)
    Public Property ModuliSpace As ModuliSpace
    Public Property RicciFlatMetric As MetricTensor
    Public Property KahlerClass As KahlerForm
    Public Property ComplexStructure As ComplexStructure
    
    ' String Theory Parameters
    Public Property StringCoupling As Double
    Public Property CompactificationScale As Double
    Public Property Supersymmetry As Integer  ' N=1,2,4,8
    
    ' Visualization Parameters
    Public Property ProjectionMatrix As Matrix4x4
    Public Property Resolution As Integer = 100
    Public Property ColorMap As Dictionary(Of Double, Color)
    
    Public Sub New()
        InitializeStandardCY()
    End Sub
    
    Private Sub InitializeStandardCY()
        ' Initialize a standard Calabi-Yau 3-fold (complex dimension 3)
        EulerCharacteristic = -200
        StringCoupling = 0.1
        CompactificationScale = 1.0E-33 ' Planck length scale
        Supersymmetry = 1
        
        ' Initialize Hodge diamond for quintic threefold
        HodgeDiamond = New Integer(3, 3) {
            {1, 0, 0, 1},
            {0, 101, 101, 0},
            {0, 101, 101, 0},
            {1, 0, 0, 1}
        }
        
        InitializeModuliSpace()
        InitializeMetric()
        InitializeColorMap()
    End Sub
    
    Private Sub InitializeModuliSpace()
        ' Moduli space = Kähler moduli × Complex structure moduli
        ModuliSpace = New ModuliSpace With {
            .KahlerDimension = HodgeDiamond(1, 1),
            .ComplexDimension = HodgeDiamond(2, 1),
            .TotalDimension = HodgeDiamond(1, 1) + HodgeDiamond(2, 1)
        }
    End Sub
    
    Private Sub InitializeMetric()
        ' Ricci-flat metric (solution to Einstein equations in vacuum)
        RicciFlatMetric = New MetricTensor(6) With {
            .IsRicciFlat = True,
            .IsKahler = True,
            .HolonomyGroup = "SU(3)"
        }
        
        ' Kähler form (symplectic structure)
        KahlerClass = New KahlerForm With {
            .IsClosed = True,
            .IsPositive = True,
            .CohomologyClass = "H^{1,1}"
        }
        
        ' Complex structure (integrable almost complex structure)
        ComplexStructure = New ComplexStructure With {
            .IsIntegrable = True,
            .NijenhuisTensor = 0,
            .Type = "Kähler"
        }
    End Sub
    
    Private Sub InitializeColorMap()
        ' Color mapping for different properties
        ColorMap = New Dictionary(Of Double, Color) From {
            {0.0, Color.Blue},
            {0.25, Color.Cyan},
            {0.5, Color.Green},
            {0.75, Color.Yellow},
            {1.0, Color.Red}
        }
    End Sub
End Class

' ==============================================
' 2. MATHEMATICAL STRUCTURES
' ==============================================

Public Structure ComplexNumber
    Public Property Real As Double
    Public Property Imaginary As Double
    
    Public Sub New(r As Double, i As Double)
        Real = r
        Imaginary = i
    End Sub
    
    Public ReadOnly Property Magnitude As Double
        Get
            Return Math.Sqrt(Real * Real + Imaginary * Imaginary)
        End Get
    End Property
    
    Public Shared Operator +(c1 As ComplexNumber, c2 As ComplexNumber) As ComplexNumber
        Return New ComplexNumber(c1.Real + c2.Real, c1.Imaginary + c2.Imaginary)
    End Operator
    
    Public Shared Operator *(c1 As ComplexNumber, c2 As ComplexNumber) As ComplexNumber
        Return New ComplexNumber(
            c1.Real * c2.Real - c1.Imaginary * c2.Imaginary,
            c1.Real * c2.Imaginary + c1.Imaginary * c2.Real
        )
    End Operator
End Structure

Public Class MetricTensor
    Public Property Dimension As Integer
    Public Property Components As Double(,)
    Public Property IsRicciFlat As Boolean
    Public Property IsKahler As Boolean
    Public Property HolonomyGroup As String
    
    Public Sub New(dim As Integer)
        Dimension = dim
        Components = New Double(dim - 1, dim - 1) {}
        InitializeAsIdentity()
    End Sub
    
    Private Sub InitializeAsIdentity()
        For i As Integer = 0 To Dimension - 1
            For j As Integer = 0 To Dimension - 1
                Components(i, j) = If(i = j, 1.0, 0.0)
            Next
        Next
    End Sub
    
    Public Function CalculateRicciTensor() As Double(,)
        ' Simplified Ricci tensor calculation
        Dim ricci As Double(,) = New Double(Dimension - 1, Dimension - 1) {}
        
        ' For Calabi-Yau, Ricci tensor should be zero
        If IsRicciFlat Then
            For i As Integer = 0 To Dimension - 1
                For j As Integer = 0 To Dimension - 1
                    ricci(i, j) = 0.0
                Next
            Next
        End If
        
        Return ricci
    End Function
End Class

Public Class KahlerForm
    Public Property IsClosed As Boolean
    Public Property IsPositive As Boolean
    Public Property CohomologyClass As String
    Public Property Volume As Double
    
    Public Function IntegrateOverCycle(cycle As HomologyCycle) As Double
        ' Integration of Kähler form over homology cycle
        Return cycle.Dimension * Math.PI
    End Function
End Class

Public Class ComplexStructure
    Public Property IsIntegrable As Boolean
    Public Property NijenhuisTensor As Double
    Public Property Type As String
    
    Public Function AlmostComplexOperator() As Double(,)
        ' J^2 = -I
        Dim J As Double(,) = New Double(5, 5) {}
        
        ' Standard complex structure for R^6 = C^3
        For i As Integer = 0 To 2
            J(2 * i, 2 * i + 1) = -1.0
            J(2 * i + 1, 2 * i) = 1.0
        Next
        
        Return J
    End Function
End Class

Public Class ModuliSpace
    Public Property KahlerDimension As Integer
    Public Property ComplexDimension As Integer
    Public Property TotalDimension As Integer
    
    Public ReadOnly Property Volume As Double
        Get
            ' Volume of moduli space (simplified)
            Return Math.Exp(TotalDimension)
        End Get
    End Property
End Class

Public Class HomologyCycle
    Public Property Dimension As Integer
    Public Property IsCalibrated As Boolean
    Public Property Volume As Double
    
    Public Sub New(dim As Integer)
        Dimension = dim
        Volume = Math.Pow(2 * Math.PI, dim)
    End Sub
End Class

' ==============================================
' 3. STRING THEORY & PHYSICS
' ==============================================

Public Class StringCompactification
    Public Property CalabiYau As CalabiYauManifold
    Public Property StringType As StringType
    Public Property GaugeGroup As String
    Public Property ChiralMatter As Integer
    Public Property GenerationCount As Integer
    
    Public Enum StringType
        TypeIIA
        TypeIIB
        HeteroticE8xE8
        HeteroticSO32
    End Enum
    
    Public Sub New(cy As CalabiYauManifold, st As StringType)
        CalabiYau = cy
        StringType = st
        CalculatePhysics()
    End Sub
    
    Private Sub CalculatePhysics()
        ' Calculate physical properties from CY geometry
        
        ' Number of generations = |χ|/2 for standard embedding
        GenerationCount = Math.Abs(CalabiYau.EulerCharacteristic) \ 2
        
        ' Determine gauge group based on string type
        Select Case StringType
            Case StringType.HeteroticE8xE8
                GaugeGroup = "E8 x E8"
                ChiralMatter = 3 * GenerationCount ' Standard Model like
                
            Case StringType.HeteroticSO32
                GaugeGroup = "SO(32)"
                ChiralMatter = 2 * GenerationCount
                
            Case StringType.TypeIIA, StringType.TypeIIB
                GaugeGroup = "U(1)^n"
                ChiralMatter = CalabiYau.HodgeDiamond(1, 1) + CalabiYau.HodgeDiamond(2, 1)
        End Select
    End Sub
    
    Public Function Calculate4DEffectiveTheory() As Effective4DTheory
        ' Compactification from 10D to 4D
        Dim theory As New Effective4DTheory
        
        ' Gravitational sector
        theory.Neutrinos = GenerationCount
        theory.Gravitinos = CalabiYau.Supersymmetry
        
        ' Gauge sector
        theory.GaugeGroup = GaugeGroup
        
        ' Matter sector
        theory.ChiralFamilies = GenerationCount
        theory.HiggsDoublets = CalabiYau.HodgeDiamond(1, 1)
        
        ' Moduli fields
        theory.KahlerModuli = CalabiYau.ModuliSpace.KahlerDimension
        theory.ComplexModuli = CalabiYau.ModuliSpace.ComplexDimension
        theory.Axions = theory.KahlerModuli + theory.ComplexModuli
        
        ' Superpotential (from flux compactification)
        theory.Superpotential = "W = ∫_CY G_3 ∧ Ω"
        
        Return theory
    End Function
End Class

Public Class Effective4DTheory
    Public Property GravitationalConstant As Double = 1.0
    Public Property CosmologicalConstant As Double = 1.0E-120
    Public Property GaugeGroup As String
    Public Property ChiralFamilies As Integer
    Public Property HiggsDoublets As Integer
    Public Property KahlerModuli As Integer
    Public Property ComplexModuli As Integer
    Public Property Axions As Integer
    Public Property Neutrinos As Integer
    Public Property Gravitinos As Integer
    Public Property Superpotential As String
End Class

' ==============================================
' 4. 3D VISUALIZATION ENGINE
' ==============================================

Public Class CalabiYauRenderer
    Inherits Control
    
    Private ReadOnly calabiYau As CalabiYauManifold
    Private ReadOnly camera As Camera
    Private ReadOnly lightSource As Vector3
    Private vertices As List(Of Vertex)
    Private edges As List(Of Edge)
    Private faces As List(Of Face)
    
    Private Class Vertex
        Public Position As Vector3
        Public Color As Color
        Public Normal As Vector3
    End Class
    
    Private Class Edge
        Public StartIndex As Integer
        Public EndIndex As Integer
        Public Color As Color
    End Class
    
    Private Class Face
        Public VertexIndices As Integer()
        Public Color As Color
        Public Normal As Vector3
    End Class
    
    Public Sub New(cy As CalabiYauManifold)
        Me.calabiYau = cy
        camera = New Camera With {
            .Position = New Vector3(0, 0, 5),
            .Target = Vector3.Zero,
            .Up = Vector3.UnitY,
            .FOV = Math.PI / 3,
            .NearPlane = 0.1,
            .FarPlane = 100
        }
        
        lightSource = New Vector3(5, 5, 5)
        
        ' Double buffering for smooth rendering
        Me.DoubleBuffered = True
        
        GenerateGeometry()
    End Sub
    
    Private Sub GenerateGeometry()
        ' Generate 3D projection of 6D Calabi-Yau manifold
        
        vertices = New List(Of Vertex)()
        edges = New List(Of Edge)()
        faces = New List(Of Face)()
        
        ' Generate points on the Calabi-Yau manifold
        Dim resolution As Integer = calabiYau.Resolution
        
        For u As Integer = 0 To resolution - 1
            For v As Integer = 0 To resolution - 1
                For w As Integer = 0 To resolution - 1
                    Dim theta1 As Double = 2 * Math.PI * u / resolution
                    Dim theta2 As Double = 2 * Math.PI * v / resolution
                    Dim theta3 As Double = 2 * Math.PI * w / resolution
                    
                    ' Map from T^6 to Calabi-Yau (simplified)
                    Dim point6D As Vector6 = MapToCalabiYau(theta1, theta2, theta3)
                    
                    ' Project to 3D for visualization
                    Dim point3D As Vector3 = ProjectTo3D(point6D)
                    
                    ' Calculate color based on curvature
                    Dim curvature As Double = CalculateGaussianCurvature(point6D)
                    Dim color As Color = GetColorFromCurvature(curvature)
                    
                    vertices.Add(New Vertex With {
                        .Position = point3D,
                        .Color = color,
                        .Normal = CalculateNormal(point6D)
                    })
                Next
            Next
        Next
        
        ' Generate edges and faces (simplified mesh)
        GenerateMesh()
    End Sub
    
    Private Function MapToCalabiYau(theta1 As Double, theta2 As Double, theta3 As Double) As Vector6
        ' Simplified mapping to Calabi-Yau manifold
        ' In reality, this would be much more complex
        
        Dim x As New Vector6
        
        ' Coordinates in 6D
        x.X = Math.Cos(theta1) + 0.1 * Math.Cos(3 * theta2)
        x.Y = Math.Sin(theta1) + 0.1 * Math.Sin(3 * theta2)
        x.Z = Math.Cos(theta2) + 0.1 * Math.Cos(3 * theta3)
        x.W = Math.Sin(theta2) + 0.1 * Math.Sin(3 * theta3)
        x.U = Math.Cos(theta3) + 0.1 * Math.Cos(3 * theta1)
        x.V = Math.Sin(theta3) + 0.1 * Math.Sin(3 * theta1)
        
        ' Apply some "Calabi-Yau-ness" (Kähler and Ricci-flat conditions)
        Dim scale As Double = 1.0 + 0.2 * Math.Sin(theta1) * Math.Sin(theta2) * Math.Sin(theta3)
        x = Vector6.Multiply(x, scale)
        
        Return x
    End Function
    
    Private Function ProjectTo3D(point6D As Vector6) As Vector3
        ' Project 6D point to 3D for visualization
        
        ' Simple orthogonal projection (could use Hopf fibration)
        Return New Vector3(
            CSng(point6D.X + 0.3 * point6D.U),
            CSng(point6D.Y + 0.3 * point6D.V),
            CSng(point6D.Z + 0.3 * point6D.W)
        )
    End Function
    
    Private Function CalculateGaussianCurvature(point6D As Vector6) As Double
        ' Simplified curvature calculation
        Dim x As Double = point6D.X
        Dim y As Double = point6D.Y
        Dim z As Double = point6D.Z
        
        ' Some oscillatory function to simulate curvature variations
        Return Math.Sin(5 * x) * Math.Cos(5 * y) * Math.Sin(5 * z)
    End Function
    
    Private Function GetColorFromCurvature(curvature As Double) As Color
        ' Map curvature to color
        Dim normalized As Double = (curvature + 1) / 2 ' Map from [-1,1] to [0,1]
        
        If calabiYau.ColorMap.ContainsKey(normalized) Then
            Return calabiYau.ColorMap(normalized)
        End If
        
        ' Interpolate between colors
        Dim keys As List(Of Double) = calabiYau.ColorMap.Keys.ToList()
        keys.Sort()
        
        For i As Integer = 0 To keys.Count - 2
            If normalized >= keys(i) AndAlso normalized <= keys(i + 1) Then
                Dim t As Double = (normalized - keys(i)) / (keys(i + 1) - keys(i))
                Return InterpolateColor(calabiYau.ColorMap(keys(i)), calabiYau.ColorMap(keys(i + 1)), t)
            End If
        Next
        
        Return Color.White
    End Function
    
    Private Function InterpolateColor(c1 As Color, c2 As Color, t As Double) As Color
        Dim r As Integer = CInt(c1.R * (1 - t) + c2.R * t)
        Dim g As Integer = CInt(c1.G * (1 - t) + c2.G * t)
        Dim b As Integer = CInt(c1.B * (1 - t) + c2.B * t)
        
        Return Color.FromArgb(255, r, g, b)
    End Function
    
    Private Function CalculateNormal(point6D As Vector6) As Vector3
        ' Calculate surface normal (simplified)
        Dim gradient As Vector6 = CalculateGradient(point6D)
        Dim normal3D As Vector3 = ProjectTo3D(gradient)
        normal3D = Vector3.Normalize(normal3D)
        
        Return normal3D
    End Function
    
    Private Function CalculateGradient(point6D As Vector6) As Vector6
        ' Numerical gradient (simplified)
        Dim epsilon As Double = 0.001
        Dim grad As New Vector6
        
        ' Calculate partial derivatives numerically
        For i As Integer = 0 To 5
            Dim pointPlus As Vector6 = point6D
            Dim pointMinus As Vector6 = point6D
            
            Select Case i
                Case 0
                    pointPlus.X += epsilon
                    pointMinus.X -= epsilon
                Case 1
                    pointPlus.Y += epsilon
                    pointMinus.Y -= epsilon
                Case 2
                    pointPlus.Z += epsilon
                    pointMinus.Z -= epsilon
                Case 3
                    pointPlus.W += epsilon
                    pointMinus.W -= epsilon
                Case 4
                    pointPlus.U += epsilon
                    pointMinus.U -= epsilon
                Case 5
                    pointPlus.V += epsilon
                    pointMinus.V -= epsilon
            End Select
            
            Dim fPlus As Double = CalculateScalarFunction(pointPlus)
            Dim fMinus As Double = CalculateScalarFunction(pointMinus)
            
            grad.SetComponent(i, (fPlus - fMinus) / (2 * epsilon))
        Next
        
        Return grad
    End Function
    
    Private Function CalculateScalarFunction(point As Vector6) As Double
        ' Some scalar function on the manifold (e.g., Kähler potential)
        Return Math.Sin(point.X) * Math.Cos(point.Y) + 
               Math.Sin(point.Z) * Math.Cos(point.W) + 
               Math.Sin(point.U) * Math.Cos(point.V)
    End Function
    
    Private Sub GenerateMesh()
        ' Generate triangle mesh from point cloud
        Dim res As Integer = calabiYau.Resolution
        
        For i As Integer = 0 To vertices.Count - res - 2
            If (i + 1) Mod res <> 0 Then ' Avoid wrapping
                ' Create two triangles for each square
                Dim v1 As Integer = i
                Dim v2 As Integer = i + 1
                Dim v3 As Integer = i + res
                Dim v4 As Integer = i + res + 1
                
                ' Triangle 1
                faces.Add(New Face With {
                    .VertexIndices = {v1, v2, v3},
                    .Color = Color.FromArgb(128, vertices(v1).Color),
                    .Normal = Vector3.Normalize(
                        (vertices(v1).Normal + vertices(v2).Normal + vertices(v3).Normal) / 3)
                })
                
                ' Triangle 2
                faces.Add(New Face With {
                    .VertexIndices = {v2, v4, v3},
                    .Color = Color.FromArgb(128, vertices(v2).Color),
                    .Normal = Vector3.Normalize(
                        (vertices(v2).Normal + vertices(v4).Normal + vertices(v3).Normal) / 3)
                })
            End If
        Next
    End Sub
    
    Protected Overrides Sub OnPaint(e As PaintEventArgs)
        MyBase.OnPaint(e)
        
        Dim g As Graphics = e.Graphics
        g.SmoothingMode = SmoothingMode.AntiAlias
        
        ' Set up transformation
        Dim viewMatrix As Matrix4x4 = Matrix4x4.CreateLookAt(
            camera.Position, camera.Target, camera.Up)
        
        Dim projectionMatrix As Matrix4x4 = Matrix4x4.CreatePerspectiveFieldOfView(
            CSng(camera.FOV), Me.Width / Me.Height, CSng(camera.NearPlane), CSng(camera.FarPlane))
        
        ' Transform and draw faces
        For Each face In faces
            Dim transformedVerts As Vector3() = New Vector3(2) {}
            
            For i As Integer = 0 To 2
                Dim vert As Vector3 = vertices(face.VertexIndices(i)).Position
                vert = Vector3.Transform(vert, viewMatrix)
                vert = Vector3.Transform(vert, projectionMatrix)
                
                ' Perspective divide and viewport transform
                vert.X = (vert.X + 1) * Me.Width / 2
                vert.Y = (-vert.Y + 1) * Me.Height / 2
                
                transformedVerts(i) = vert
            Next
            
            ' Calculate lighting
            Dim lightDir As Vector3 = Vector3.Normalize(lightSource - vertices(face.VertexIndices(0)).Position)
            Dim intensity As Single = Math.Max(0, Vector3.Dot(face.Normal, lightDir))
            
            Dim litColor As Color = Color.FromArgb(
                face.Color.A,
                CInt(face.Color.R * intensity),
                CInt(face.Color.G * intensity),
                CInt(face.Color.B * intensity))
            
            ' Draw triangle
            Using brush As New SolidBrush(litColor)
                Dim points As PointF() = {
                    New PointF(transformedVerts(0).X, transformedVerts(0).Y),
                    New PointF(transformedVerts(1).X, transformedVerts(1).Y),
                    New PointF(transformedVerts(2).X, transformedVerts(2).Y)
                }
                g.FillPolygon(brush, points)
            End Using
        Next
        
        ' Draw info overlay
        DrawInfoOverlay(g)
    End Sub
    
    Private Sub DrawInfoOverlay(g As Graphics)
        Dim font As New Font("Arial", 10)
        Dim brush As New SolidBrush(Color.White)
        
        Dim info As String = String.Format(
            "Calabi-Yau 3-fold (6D Compactification)" & vbCrLf &
            "Euler Characteristic: {0}" & vbCrLf &
            "Hodge Numbers: h^{{1,1}}={1}, h^{{2,1}}={2}" & vbCrLf &
            "Moduli Space: {3}D (Kähler) + {4}D (Complex)" & vbCrLf &
            "Supersymmetry: N={5}",
            calabiYau.EulerCharacteristic,
            calabiYau.HodgeDiamond(1, 1),
            calabiYau.HodgeDiamond(2, 1),
            calabiYau.ModuliSpace.KahlerDimension,
            calabiYau.ModuliSpace.ComplexDimension,
            calabiYau.Supersymmetry)
        
        g.DrawString(info, font, brush, New PointF(10, 10))
    End Sub
End Class

Public Class Camera
    Public Property Position As Vector3
    Public Property Target As Vector3
    Public Property Up As Vector3
    Public Property FOV As Double
    Public Property NearPlane As Double
    Public Property FarPlane As Double
End Class

Public Structure Vector6
    Public X, Y, Z, W, U, V As Double
    
    Public Shared Operator *(v As Vector6, scalar As Double) As Vector6
        Return New Vector6 With {
            .X = v.X * scalar,
            .Y = v.Y * scalar,
            .Z = v.Z * scalar,
            .W = v.W * scalar,
            .U = v.U * scalar,
            .V = v.V * scalar
        }
    End Operator
    
    Public Sub SetComponent(index As Integer, value As Double)
        Select Case index
            Case 0 : X = value
            Case 1 : Y = value
            Case 2 : Z = value
            Case 3 : W = value
            Case 4 : U = value
            Case 5 : V = value
        End Select
    End Sub
End Structure

' ==============================================
' 5. MAIN FORM & USER INTERFACE
' ==============================================

Public Class CalabiYauExplorerForm
    Inherits Form
    
    Private ReadOnly renderer As CalabiYauRenderer
    Private ReadOnly calabiYau As CalabiYauManifold
    Private ReadOnly stringTheory As StringCompactification
    Private ReadOnly timer As Timer
    
    Private panel3D As Panel
    Private lstProperties As ListBox
    Private btnRotate As Button
    Private btnChangeView As Button
    Private btnExport As Button
    Private lblInfo As Label
    
    Public Sub New()
        Me.Text = "Calabi-Yau Manifold Explorer - 6D Compactification"
        Me.Size = New Size(1200, 800)
        Me.StartPosition = FormStartPosition.CenterScreen
        
        ' Initialize Calabi-Yau manifold
        calabiYau = New CalabiYauManifold()
        
        ' Initialize string theory compactification
        stringTheory = New StringCompactification(calabiYau, StringCompactification.StringType.HeteroticE8xE8)
        
        ' Setup UI
        InitializeComponents()
        
        ' Initialize renderer
        renderer = New CalabiYauRenderer(calabiYau)
        panel3D.Controls.Add(renderer)
        renderer.Dock = DockStyle.Fill
        
        ' Setup animation timer
        timer = New Timer With {.Interval = 16} ' ~60 FPS
        AddHandler timer.Tick, AddressOf Timer_Tick
        timer.Start()
        
        UpdatePropertiesList()
    End Sub
    
    Private Sub InitializeComponents()
        ' Main split container
        Dim splitContainer As New SplitContainer With {
            .Dock = DockStyle.Fill,
            .SplitterDistance = 800
        }
        
        ' 3D visualization panel
        panel3D = New Panel With {
            .Dock = DockStyle.Fill,
            .BackColor = Color.Black
        }
        
        ' Properties panel
        Dim panelProps As New Panel With {
            .Dock = DockStyle.Fill,
            .BackColor = Color.FromArgb(30, 30, 30)
        }
        
        ' Properties list
        lstProperties = New ListBox With {
            .Dock = DockStyle.Top,
            .Height = 400,
            .BackColor = Color.FromArgb(40, 40, 40),
            .ForeColor = Color.White,
            .Font = New Font("Consolas", 10)
        }
        
        ' Information label
        lblInfo = New Label With {
            .Dock = DockStyle.Top,
            .Height = 200,
            .ForeColor = Color.LightGreen,
            .BackColor = Color.FromArgb(20, 20, 20),
            .Font = New Font("Arial", 9),
            .Padding = New Padding(10)
        }
        
        ' Buttons panel
        Dim panelButtons As New FlowLayoutPanel With {
            .Dock = DockStyle.Bottom,
            .Height = 50,
            .FlowDirection = FlowDirection.LeftToRight
        }
        
        btnRotate = New Button With {
            .Text = "Auto-Rotate",
            .Size = New Size(100, 30)
        }
        AddHandler btnRotate.Click, AddressOf BtnRotate_Click
        
        btnChangeView = New Button With {
            .Text = "Change View",
            .Size = New Size(100, 30)
        }
        AddHandler btnChangeView.Click, AddressOf BtnChangeView_Click
        
        btnExport = New Button With {
            .Text = "Export Data",
            .Size = New Size(100, 30)
        }
        AddHandler btnExport.Click, AddressOf BtnExport_Click
        
        ' Add controls
        panelButtons.Controls.AddRange({btnRotate, btnChangeView, btnExport})
        
        panelProps.Controls.AddRange({lstProperties, lblInfo, panelButtons})
        
        splitContainer.Panel1.Controls.Add(panel3D)
        splitContainer.Panel2.Controls.Add(panelProps)
        
        Me.Controls.Add(splitContainer)
    End Sub
    
    Private Sub UpdatePropertiesList()
        lstProperties.Items.Clear()
        
        ' Add Calabi-Yau properties
        lstProperties.Items.Add("=== CALABI-YAU MANIFOLD ===")
        lstProperties.Items.Add($"Dimension: {calabiYau.Dimension}D")
        lstProperties.Items.Add($"Euler Characteristic: {calabiYau.EulerCharacteristic}")
        lstProperties.Items.Add($"Hodge Diamond:")
        lstProperties.Items.Add($"  h^{{0,0}} = {calabiYau.HodgeDiamond(0, 0)}")
        lstProperties.Items.Add($"  h^{{1,1}} = {calabiYau.HodgeDiamond(1, 1)}")
        lstProperties.Items.Add($"  h^{{2,1}} = {calabiYau.HodgeDiamond(2, 1)}")
        lstProperties.Items.Add($"  h^{{3,0}} = {calabiYau.HodgeDiamond(3, 0)}")
        
        lstProperties.Items.Add("")
        lstProperties.Items.Add("=== MODULI SPACE ===")
        lstProperties.Items.Add($"Kähler Moduli: {calabiYau.ModuliSpace.KahlerDimension}")
        lstProperties.Items.Add($"Complex Moduli: {calabiYau.ModuliSpace.ComplexDimension}")
        lstProperties.Items.Add($"Total Dimension: {calabiYau.ModuliSpace.TotalDimension}")
        
        ' Add string theory properties
        Dim theory As Effective4DTheory = stringTheory.Calculate4DEffectiveTheory()
        
        lstProperties.Items.Add("")
        lstProperties.Items.Add("=== STRING COMPACTIFICATION ===")
        lstProperties.Items.Add($"String Type: {stringTheory.StringType}")
        lstProperties.Items.Add($"Gauge Group: {theory.GaugeGroup}")
        lstProperties.Items.Add($"Chiral Families: {theory.ChiralFamilies}")
        lstProperties.Items.Add($"Generations: {stringTheory.GenerationCount}")
        lstProperties.Items.Add($"Supersymmetry: N={calabiYau.Supersymmetry}")
        
        ' Update info label
        lblInfo.Text = $"Compactification Scale: {calabiYau.CompactificationScale:E3} m" & vbCrLf &
                       $"String Coupling: {calabiYau.StringCoupling}" & vbCrLf &
                       $"Axions: {theory.Axions}" & vbCrLf &
                       $"Superpotential: {theory.Superpotential}"
    End Sub
    
    Private rotationAngle As Double = 0
    
    Private Sub Timer_Tick(sender As Object, e As EventArgs)
        ' Auto-rotate the visualization
        rotationAngle += 0.01
        renderer.Invalidate() ' Trigger repaint
    End Sub
    
    Private Sub BtnRotate_Click(sender As Object, e As EventArgs)
        ' Toggle auto-rotation
        If timer.Enabled Then
            timer.Stop()
            btnRotate.Text = "Start Rotation"
        Else
            timer.Start()
            btnRotate.Text = "Stop Rotation"
        End If
    End Sub
    
    Private Sub BtnChangeView_Click(sender As Object, e As EventArgs)
        ' Cycle through different views
        Static viewIndex As Integer = 0
        viewIndex = (viewIndex + 1) Mod 4
        
        Select Case viewIndex
            Case 0 ' Standard view
                renderer.BackColor = Color.Black
            Case 1 ' Wireframe view
                renderer.BackColor = Color.DarkBlue
            Case 2 ' Curvature map
                renderer.BackColor = Color.DarkGray
            Case 3 ' Topological view
                renderer.BackColor = Color.DarkGreen
        End Select
    End Sub
    
    Private Sub BtnExport_Click(sender As Object, e As EventArgs)
        ' Export CY data
        Using dialog As New SaveFileDialog
            dialog.Filter = "Text Files|*.txt|XML Files|*.xml|JSON Files|*.json"
            dialog.Title = "Export Calabi-Yau Data"
            
            If dialog.ShowDialog() = DialogResult.OK Then
                ExportData(dialog.FileName)
                MessageBox.Show($"Data exported to {dialog.FileName}", "Export Complete", 
                              MessageBoxButtons.OK, MessageBoxIcon.Information)
            End If
        End Using
    End Sub
    
    Private Sub ExportData(filename As String)
        Dim data As String = $"Calabi-Yau Manifold Data Export" & vbCrLf &
                           $"==============================" & vbCrLf &
                           $"Timestamp: {DateTime.Now}" & vbCrLf &
                           $"Euler Characteristic: {calabiYau.EulerCharacteristic}" & vbCrLf &
                           $"Hodge Numbers:" & vbCrLf &
                           $"  h^(1,1) = {calabiYau.HodgeDiamond(1, 1)}" & vbCrLf &
                           $"  h^(2,1) = {calabiYau.HodgeDiamond(2, 1)}" & vbCrLf &
                           $"Moduli Space Dimensions: {calabiYau.ModuliSpace.TotalDimension}" & vbCrLf &
                           $"String Coupling: {calabiYau.StringCoupling}" & vbCrLf &
                           $"Compactification Scale: {calabiYau.CompactificationScale} m"
        
        System.IO.File.WriteAllText(filename, data)
    End Sub
    
    Protected Overrides Sub OnKeyDown(e As KeyEventArgs)
        MyBase.OnKeyDown(e)
        
        ' Keyboard controls for camera
        Select Case e.KeyCode
            Case Keys.Left
                rotationAngle -= 0.1
            Case Keys.Right
                rotationAngle += 0.1
            Case Keys.Up
                ' Zoom in
            Case Keys.Down
                ' Zoom out
            Case Keys.R
                ' Reset view
                rotationAngle = 0
        End Select
        
        renderer.Invalidate()
    End Sub
End Class

' ==============================================
' 6. PROGRAM ENTRY POINT
' ==============================================

Module Program
    <STAThread>
    Sub Main()
        Application.EnableVisualStyles()
        Application.SetCompatibleTextRenderingDefault(False)
        Application.Run(New CalabiYauExplorerForm())
    End Sub
End Module