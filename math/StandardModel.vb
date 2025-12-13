Imports System.Collections.Generic

Public Enum Particle
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
    Public Property Name As String
    Public Property Type As String
    Public Property Charge As String
    Public Property Spin As String
    Public Property Mass As String

    Public Sub New(name As String, type As String, charge As String, spin As String, mass As String)
        Me.Name = name
        Me.Type = type
        Me.Charge = charge
        Me.Spin = spin
        Me.Mass = mass
    End Sub
End Class

Module ParticleData
    Public ReadOnly InfoMap As New Dictionary(Of Particle, Info) From {
        {Particle.UpQuark, New Info("Up Quark", "Quark", "+2/3", "1/2", "2.2 MeV/c²")},
        {Particle.CharmQuark, New Info("Charm Quark", "Quark", "+2/3", "1/2", "1.28 GeV/c²")},
        {Particle.TopQuark, New Info("Top Quark", "Quark", "+2/3", "1/2", "173.1 GeV/c²")},
        {Particle.DownQuark, New Info("Down Quark", "Quark", "-1/3", "1/2", "4.7 MeV/c²")},
        {Particle.StrangeQuark, New Info("Strange Quark", "Quark", "-1/3", "1/2", "96 MeV/c²")},
        {Particle.BottomQuark, New Info("Bottom Quark", "Quark", "-1/3", "1/2", "4.18 GeV/c²")},
        {Particle.Electron, New Info("Electron", "Lepton", "-1", "1/2", "0.511 MeV/c²")},
        {Particle.MuonLepton, New Info("Muon", "Lepton", "-1", "1/2", "105.7 MeV/c²")},
        {Particle.TauLepton, New Info("Tau", "Lepton", "-1", "1/2", "1.777 GeV/c²")},
        {Particle.ElectronNeutrino, New Info("Electron Neutrino", "Lepton", "0", "1/2", "<2.2 eV/c²")},
        {Particle.MuonNeutrino, New Info("Muon Neutrino", "Lepton", "0", "1/2", "<0.17 MeV/c²")},
        {Particle.TauNeutrino, New Info("Tau Neutrino", "Lepton", "0", "1/2", "<18.2 MeV/c²")},
        {Particle.Photon, New Info("Photon", "Gauge Boson", "0", "1", "0")},
        {Particle.GluonBoson, New Info("Gluon", "Gauge Boson", "0", "1", "0")},
        {Particle.ZBoson, New Info("Z Boson", "Gauge Boson", "0", "1", "91.2 GeV/c²")},
        {Particle.WPlusBoson, New Info("W+ Boson", "Gauge Boson", "+1", "1", "80.4 GeV/c²")},
        {Particle.WMinusBoson, New Info("W- Boson", "Gauge Boson", "-1", "1", "80.4 GeV/c²")},
        {Particle.HiggsBoson, New Info("Higgs Boson", "Scalar Boson", "0", "0", "125.1 GeV/c²")}
    }
End Module