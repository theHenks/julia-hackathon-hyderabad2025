using UnROOT
using LorentzVectorHEP
using CairoMakie

# Get the CMS OpenData file from CERN EOS
#file = "root://eospublic.cern.ch//eos/opendata/cms/Run2016H/MuonEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/130000/0E3D51CB-3B75-974F-B868-0E2ABA272073.root"
# Alternatively, use a local file path:
file = "/data/cms/0E3D51CB-3B75-974F-B868-0E2ABA272073.root"

# Open the ROOT file and access the Events tree selecting all branches matching "Muon.*"
tfile = ROOTFile(file)         
events = LazyTree(tfile, "Events",  Regex("Muon.*"));

# Loop over events and select dimuon candidates with specific criteria 
masses = Float32[]
for evt in events
    evt.nMuon == 2 || continue                  # exactly two muons
    sum(evt.Muon_charge) == 0 || continue       # opposite sign
    all(>(10.), evt.Muon_pt) || continue        # pT > 10 GeV
    muons = LorentzVectorCyl.(evt.Muon_pt, evt.Muon_eta, evt.Muon_phi, evt.Muon_mass)
    mass = fast_mass(muons[1], muons[2]) 
    (60.0 < mass < 120.0) || continue           # Z mass window
    push!(masses, mass)
end

# plot the dimuon mass spectrum
fig = Figure(size=(600, 400))
hist(fig[1,1], masses, bins=100, axis = (xlabel="Dimuon Mass [GeV]", ylabel="Events", title="Dimuon Invariant Mass Spectrum"))
fig
