Imports System
Imports System.Collections.Generic
Imports System.Text

''' Represents a solid body (convex polyhedron) based on Euler's Polyhedron Equation V - E + F = 2.vertices, edges, and faces.
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

    ''' Euler Characteristic constant for convex polyhedra
    Public Const EULER_CHARACTERISTIC As Integer = 2

    ''' 3D Point structure for vertex representation
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
    
    ''' Edge structure connecting two vertices
    
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

    ''' Face structure defined by vertex indices
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
    Public Sub New(name, vertices, edges, faces, Optional faceType As String = "Mixed")
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

    ''' Euler's Polyhedron Formula V - E + F = 2
    Public Function ValidateEulerFormula() As Boolean
        Return (_vertices - _edges + _faces) = EULER_CHARACTERISTIC
    End Function
    
    ''' Calculates 
    Public ReadOnly Property EulerCharacteristic As Integer
        Get
            Return _vertices - _edges + _faces
        End Get
    End Property
    
    ''' Verifies if the polyhedron satisfies necessary conditions
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
    
    ''' Calculates the average number of edges per face
    Public ReadOnly Property EdgesPerFace As Double
        Get
            If _faces = 0 Then Return 0
            Return (2.0 * _edges) / _faces
        End Get
    End Property

    ''' Calculates the average number of edges meeting at each vertex           
    Public ReadOnly Property EdgesPerVertex As Double
        Get
            If _vertices = 0 Then Return 0
            Return (2.0 * _edges) / _vertices
        End Get
    End Property
    
    ''' Generates a report about the solid body
    Public Function GetEulerReport() As String
        Dim sb As New StringBuilder()
        
        Return sb.ToString()
    End Function
    
    ''' Adds a vertex to the solid body
    Sub AddVertex(point As Point3D)
        _vertexList.Add(point)
    End Sub
    
    ''' Adds an edge to the solid body
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
                                
''' Factory class for creating common Platonic solids
Public Class PlatonicSolidFactory
    ''' Creates a tetrahedron (4 faces, all triangles)
    Shared Function CreateTetrahedron() As SolidBody
        Dim tetra = New SolidBody("Tetrahedron", 4, 6, 4, "Triangle")
        ' Note: In a full implementation, you would add actual vertices, edges, and faces here
        Return tetra
    End Function
    ''' Creates a cube/hexahedron (6 faces, all squares)
    Shared Function Cube() As SolidBody
        Return New SolidBody("Cube", 8, 12, 6, "Square")
    End Function
    ''' Creates an octahedron (8 faces, all triangles)
    Shared Function Octahedron() As SolidBody
        Return New SolidBody("Octahedron", 6, 12, 8, "Triangle")
    End Function
    ''' Creates a dodecahedron (12 faces, all pentagons)
    Shared Function Dodecahedron() As SolidBody
        Return New SolidBody("Dodecahedron", 20, 30, 12, "Pentagon")
    End Function
    ''' Creates an icosahedron (20 faces, all triangles)
    Shared Function Icosahedron() As SolidBody
        Return New SolidBody("Icosahedron", 12, 30, 20, "Triangle")
    End Function
    ''' Returns all five Platonic solids
    Shared Function GetAllPlatonicSolids() As List(Of SolidBody)
        Return New List(Of SolidBody) From {
            Tetrahedron(),
            Cube(),
            Octahedron(),
            Dodecahedron(),
            Icosahedron()
        }
    End Function
End Class
