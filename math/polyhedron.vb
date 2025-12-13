Imports System
Imports System.Collections.Generic
Imports System.Text

''' <summary>
''' Represents a solid body (convex polyhedron) based on Euler's Polyhedron Equation V - E + F = 2.
''' Models a three-dimensional shape with vertices, edges, and faces.
''' </summary>
Public Class SolidBody
    ' Core properties representing Euler's equation components
    Private _vertices As Integer
    Private _edges As Integer
    Private _faces As Integer
    
    ' Additional properties
    Private _name As String
    Private _faceType As String
    Private _isValid As Boolean
    
    ' Collections for detailed modeling
    Private _vertexList As List(Of Point3D)
    Private _edgeList As List(Of Edge)
    Private _faceList As List(Of Face)
    
    ''' <summary>
    ''' Euler Characteristic constant for convex polyhedra
    ''' </summary>
    Public Const EULER_CHARACTERISTIC As Integer = 2

    ''' <summary>
    ''' 3D Point structure for vertex representation
    ''' </summary>
    Public Structure Point3D
        Public X As Double
        Public Y As Double
        Public Z As Double
        
        Public Sub New(x As Double, y As Double, z As Double)
            Me.X = x
            Me.Y = y
            Me.Z = z
        End Sub
        
        Public Overrides Function ToString() As String
            Return $"({X:F2}, {Y:F2}, {Z:F2})"
        End Function
    End Structure

    ''' <summary>
    ''' Edge structure connecting two vertices
    ''' </summary>
    Public Class Edge
        Public Property StartVertex As Integer
        Public Property EndVertex As Integer
        Public Property Length As Double
        
        Public Sub New(start As Integer, [end] As Integer, length As Double)
            Me.StartVertex = start
            Me.EndVertex = [end]
            Me.Length = length
        End Sub
        
        Public Overrides Function ToString() As String
            Return $"Edge V{StartVertex} → V{EndVertex} (Length: {Length:F2})"
        End Function
    End Class

    ''' <summary>
    ''' Face structure defined by vertex indices
    ''' </summary>
    Public Class Face
        Public Property VertexIndices As List(Of Integer)
        Public Property FaceType As String
        Public Property Area As Double
        
        Public Sub New(vertices As List(Of Integer), Optional faceType As String = "Polygon")
            Me.VertexIndices = vertices
            Me.FaceType = faceType
            Me.Area = 0
        End Sub
        
        Public ReadOnly Property VertexCount As Integer
            Get
                Return VertexIndices.Count
            End Get
        End Property
        
        Public Overrides Function ToString() As String
            Return $"{FaceType} face with {VertexCount} vertices"
        End Function
    End Class

    ''' <summary>
    ''' Creates a new SolidBody with specified parameters
    ''' </summary>
    ''' <param name="name">Name of the polyhedron</param>
    ''' <param name="vertices">Number of vertices (V)</param>
    ''' <param name="edges">Number of edges (E)</param>
    ''' <param name="faces">Number of faces (F)</param>
    ''' <param name="faceType">Type of faces (e.g., "Triangle", "Square")</param>
    Public Sub New(name As String, vertices As Integer, edges As Integer, faces As Integer, Optional faceType As String = "Mixed")
        _name = name
        _vertices = vertices
        _edges = edges
        _faces = faces
        _faceType = faceType
        
        ' Validate against Euler's formula
        _isValid = ValidateEulerFormula()
        
        ' Initialize collections
        _vertexList = New List(Of Point3D)()
        _edgeList = New List(Of Edge)()
        _faceList = New List(Of Face)()
    End Sub
    
    ''' <summary>
    ''' Validates the solid body using Euler's Polyhedron Formula V - E + F = 2
    ''' </summary>
    Public Function ValidateEulerFormula() As Boolean
        Return (_vertices - _edges + _faces) = EULER_CHARACTERISTIC
    End Function
    
    ''' <summary>
    ''' Calculates the Euler Characteristic
    ''' </summary>
    Public ReadOnly Property EulerCharacteristic As Integer
        Get
            Return _vertices - _edges + _faces
        End Get
    End Property
    
    ''' <summary>
    ''' Verifies if the polyhedron satisfies necessary conditions
    ''' </summary>
    Public Function IsValidPolyhedron() As Boolean
        ' Check Euler's formula
        If Not ValidateEulerFormula() Then Return False
        
        ' Additional polyhedron constraints:
        ' 1. All counts must be positive
        If _vertices <= 0 OrElse _edges <= 0 OrElse _faces <= 0 Then Return False
        
        ' 2. Edge-Vertex relationship: Each edge connects 2 vertices
        If _edges > (_vertices * (_vertices - 1)) / 2 Then Return False
        
        ' 3. Face-Vertex relationship: Minimum 4 faces for tetrahedron
        If _faces < 4 Then Return False
        
        Return True
    End Function
    
    ''' <summary>
    ''' Calculates the average number of edges per face
    ''' </summary>
    Public ReadOnly Property EdgesPerFace As Double
        Get
            If _faces = 0 Then Return 0
            Return (2.0 * _edges) / _faces
        End Get
    End Property
    
    ''' <summary>
    ''' Calculates the average number of edges meeting at each vertex
    ''' </summary>
    Public ReadOnly Property EdgesPerVertex As Double
        Get
            If _vertices = 0 Then Return 0
            Return (2.0 * _edges) / _vertices
        End Get
    End Property
    
    ''' <summary>
    ''' Generates a report about the solid body
    ''' </summary>
    Public Function GetEulerReport() As String
        Dim sb As New StringBuilder()
        
        sb.AppendLine($"Solid Body: {_name}")
        sb.AppendLine($"{"=".PadRight(40, "=")}")
        sb.AppendLine($"Vertices (V): {_vertices}")
        sb.AppendLine($"Edges (E): {_edges}")
        sb.AppendLine($"Faces (F): {_faces}")
        sb.AppendLine()
        sb.AppendLine($"Euler Calculation: V - E + F = {_vertices} - {_edges} + {_faces}")
        sb.AppendLine($"Result: {EulerCharacteristic}")
        sb.AppendLine($"Valid according to Euler: {ValidateEulerFormula()}")
        sb.AppendLine()
        sb.AppendLine($"Topological Type: {(If(ValidateEulerFormula(), "Sphere (χ=2)", "Non-spherical topology"))}")
        sb.AppendLine($"Average edges per face: {EdgesPerFace:F2}")
        sb.AppendLine($"Average edges per vertex: {EdgesPerVertex:F2}")
        
        Return sb.ToString()
    End Function
    
    ' Properties with getters
    Public ReadOnly Property Name As String
        Get
            Return _name
        End Get
    End Property
    
    Public ReadOnly Property Vertices As Integer
        Get
            Return _vertices
        End Get
    End Property
    
    Public ReadOnly Property Edges As Integer
        Get
            Return _edges
        End Get
    End Property
    
    Public ReadOnly Property Faces As Integer
        Get
            Return _faces
        End Get
    End Property
    
    Public ReadOnly Property FaceType As String
        Get
            Return _faceType
        End Get
    End Property
    
    Public ReadOnly Property IsValid As Boolean
        Get
            Return _isValid
        End Get
    End Property
    
    ' Collection accessors
    Public ReadOnly Property VertexPoints As List(Of Point3D)
        Get
            Return _vertexList
        End Get
    End Property
    
    Public ReadOnly Property EdgesList As List(Of Edge)
        Get
            Return _edgeList
        End Get
    End Property
    
    Public ReadOnly Property FacesList As List(Of Face)
        Get
            Return _faceList
        End Get
    End Property
    
    ''' <summary>
    ''' Adds a vertex to the solid body
    ''' </summary>
    Public Sub AddVertex(point As Point3D)
        _vertexList.Add(point)
    End Sub
    
    ''' <summary>
    ''' Adds an edge to the solid body
    ''' </summary>
    Public Sub AddEdge(edge As Edge)
        _edgeList.Add(edge)
    End Sub
    
    ''' <summary>
    ''' Adds a face to the solid body
    ''' </summary>
    Public Sub AddFace(face As Face)
        _faceList.Add(face)
    End Sub
End Class

''' <summary>
''' Factory class for creating common Platonic solids
''' </summary>
Public Class PlatonicSolidFactory
    ''' <summary>
    ''' Creates a tetrahedron (4 faces, all triangles)
    ''' </summary>
    Public Shared Function CreateTetrahedron() As SolidBody
        Dim tetra = New SolidBody("Tetrahedron", 4, 6, 4, "Triangle")
        ' Note: In a full implementation, you would add actual vertices, edges, and faces here
        Return tetra
    End Function
    
    ''' <summary>
    ''' Creates a cube/hexahedron (6 faces, all squares)
    ''' </summary>
    Public Shared Function CreateCube() As SolidBody
        Return New SolidBody("Cube", 8, 12, 6, "Square")
    End Function
    
    ''' <summary>
    ''' Creates an octahedron (8 faces, all triangles)
    ''' </summary>
    Public Shared Function CreateOctahedron() As SolidBody
        Return New SolidBody("Octahedron", 6, 12, 8, "Triangle")
    End Function
    
    ''' <summary>
    ''' Creates a dodecahedron (12 faces, all pentagons)
    ''' </summary>
    Public Shared Function CreateDodecahedron() As SolidBody
        Return New SolidBody("Dodecahedron", 20, 30, 12, "Pentagon")
    End Function
    
    ''' <summary>
    ''' Creates an icosahedron (20 faces, all triangles)
    ''' </summary>
    Public Shared Function CreateIcosahedron() As SolidBody
        Return New SolidBody("Icosahedron", 12, 30, 20, "Triangle")
    End Function
    
    ''' <summary>
    ''' Returns all five Platonic solids
    ''' </summary>
    Public Shared Function GetAllPlatonicSolids() As List(Of SolidBody)
        Return New List(Of SolidBody) From {
            CreateTetrahedron(),
            CreateCube(),
            CreateOctahedron(),
            CreateDodecahedron(),
            CreateIcosahedron()
        }
    End Function
End Class

''' <summary>
''' Main module demonstrating usage of the SolidBody class
''' </summary>
Module MainModule
    Sub Main()
        Console.WriteLine("SOLID BODY ANALYSIS USING EULER'S FORMULA V - E + F = 2")
        Console.WriteLine("=".PadRight(60, "="))
        Console.WriteLine()
        
        ' Example 1: Create and analyze individual solids
        Dim cube As SolidBody = PlatonicSolidFactory.CreateCube()
        Console.WriteLine(cube.GetEulerReport())
        
        Console.WriteLine("-".PadRight(40, "-"))
        
        Dim tetra As SolidBody = PlatonicSolidFactory.CreateTetrahedron()
        Console.WriteLine(tetra.GetEulerReport())
        
        Console.WriteLine()
        Console.WriteLine("ALL PLATONIC SOLIDS:")
        Console.WriteLine("-".PadRight(40, "-"))
        
        ' Example 2: Analyze all Platonic solids
        Dim allSolids = PlatonicSolidFactory.GetAllPlatonicSolids()
        
        For Each solid In allSolids
            Console.WriteLine($"{solid.Name.PadRight(15)}: V={solid.Vertices.ToString().PadLeft(2)}, " &
                             $"E={solid.Edges.ToString().PadLeft(2)}, F={solid.Faces.ToString().PadLeft(2)} | " &
                             $"V-E+F={solid.EulerCharacteristic} | Valid: {solid.IsValid}")
        Next
        
        Console.WriteLine()
        Console.WriteLine("VERIFICATION OF EULER'S FORMULA:")
        Console.WriteLine("-".PadRight(40, "-"))
        
        ' Example 3: Test with invalid polyhedron
        Dim invalidSolid As New SolidBody("Invalid Shape", 6, 9, 4, "Triangle")
        Console.WriteLine(invalidSolid.GetEulerReport())
        Console.WriteLine($"Is valid polyhedron: {invalidSolid.IsValidPolyhedron()}")
        
        Console.WriteLine()
        Console.WriteLine("ADDITIONAL PROPERTIES:")
        Console.WriteLine("-".PadRight(40, "-"))
        
        ' Example 4: Calculate derived properties
        Dim octahedron = PlatonicSolidFactory.CreateOctahedron()
        Console.WriteLine($"Octahedron analysis:")
        Console.WriteLine($"  Average edges per face: {octahedron.EdgesPerFace:F2}")
        Console.WriteLine($"  Average edges per vertex: {octahedron.EdgesPerVertex:F2}")
        Console.WriteLine($"  Face type: {octahedron.FaceType}")
        
        Console.WriteLine()
        Console.WriteLine("Press any key to exit...")
        Console.ReadKey()
    End Sub
End Module