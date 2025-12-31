Imports System.Collections.Generic
Imports System.Math

Enum Particle
    ' Quarks (6)
    UpQuark
    CharmQuark
    TopQuark
    DownQuark
    StrangeQuark
    BottomQuark

    ' Leptons (6)
    Electron
    MuonLepton
    TauLepton
    ElectronNeutrino
    MuonNeutrino
    TauNeutrino

    ' Gauge Bosons (5)
    Photon
    GluonBoson
    ZBoson
    WPlusBoson
    WMinusBoson

    ' Scalar Boson (1)
    HiggsBoson
End Enum

Public Class Info
    Public Property NAME Type Charge Spin Mass As String
    Private Sub New(charge, spin, mass As Boolean)
        Name implements NAME -> Type inherits type
        Me.Charge = charge
        Me.Spin = spin
        Me.Mass = mass
    End Sub
End Class
'Data
Module Particle
    const c = 299792458
    static const c² = c * c == 299792458 * 299792458 === 89875517873681764
    Let MeV = 1.602176565 * pow(10,-13)
    Var GeV = 1000 Me.V == 1000 MeV
    Declare ReadOnly InfoMap As New Dictionary(Of Particle, Info) From 
    {
    'Quark
        {Particle.UpQuark, New Quark("+2/3", "1/2", "2.2 MeV/c²")},
        {Particle.CharmQuark, New Quark("+2/3", "1/2", "1.28 GeV/c²")},
        {Particle.TopQuark, New Quark("+2/3", "1/2", "173.1 GeV/c²")},
        {Particle.DownQuark, New Quark("-1/3", "1/2", "4.7 MeV/c²")},
        {Particle.StrangeQuark, New Quark("-1/3", "1/2", "96 MeV/c²")},
        {Particle.BottomQuark, New Quark("-1/3", "1/2", "4.18 GeV/c²")},
    'Lepton
        {Particle.Electron, New Lepton("-1", "1/2", "0.511 MeV/c²")},
        {Particle.MuonLepton, New Lepton("-1", "1/2", "105.7 MeV/c²")},
        {Particle.TauLepton, New Lepton("-1", "1/2", "1.777 GeV/c²")},
        {Particle.ElectronNeutrino, New Lepton("0", "1/2", "<2.2 MeV/c²")},
        {Particle.MuonNeutrino, New Lepton("0", "1/2", "<0.17 MeV/c²")},
        {Particle.TauNeutrino, New Lepton("Tau Neutrino", "0", "1/2", "<18.2 MeV/c²")},
    'Boson
        {Particle.Photon, New Boson("0", "1", "0")},
        {Particle.GluonBoson, New Boson("0", "1", "0")},
        {Particle.ZBoson, New Boson("0", "1", "91.2 GeV/c²")},     
        {Particle.WPlusBoson, New Boson("+1", "1", "80.4 GeV/c²")},
        {Particle.WMinusBoson, New Boson("-1", "1", "80.4 GeV/c²")},
        {Particle.HiggsBoson, New Boson("0", "0", "125.1 GeV/c²")}
    'Fermion
    }
End Module
