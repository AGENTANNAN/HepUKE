# Paper: 2206.13674v1
# Title: Measurement of the total and leptonic decay widths of the J/psi resonance
#        with an energy scan method at BESIII
# Energy: 3.050-3.120 GeV (15 energy points, J/psi scan)
# Processes: e+e- -> e+e- (Bhabha) and e+e- -> mu+mu- (dimuon)
# QED cross-section measurement, no hadronic decays
# ConExc continuum analysis

### Dataset preparation ###
# 15 energy points around J/psi: 3049.642 - 3119.878 MeV
# Using inclusive MC for background estimation
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> e+e- (Bhabha scattering)
decay_card_ee = <<~DECAYCARD
    Decay e+ e-
    1.000  e+  e-                                 CONEXC;
    Enddecay
End
DECAYCARD

# Decay card: e+e- -> mu+mu- (dimuon production)
decay_card_mumu = <<~DECAYCARD
    Decay e+ e-
    1.000  mu+  mu-                               CONEXC;
    Enddecay
End
DECAYCARD

# Exclusive MC for both processes
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "ee_to_ee_Bhabha"
  config.events         = 500000
  config.decay_card     = decay_card_ee
  config.cross_section  = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "ee_to_mumu"
  config.events         = 500000
  config.decay_card     = decay_card_mumu
  config.cross_section  = :default
end

### Event selection: e+e- -> e+e- channel ###
alg_ee = Algorithm.new("Jpsi_scan_ee")
alg_ee.set_header(["JpsiScanEEAlg/JpsiScanEE.h"])

event_selection_ee = Selection.new

# Two oppositely charged tracks, barrel region
event_selection_ee.select_track {
                   cos_theta 0.8
                   Vz   10.0
                   Vr   1.0
                   nTot "==2"
                 }
                 .assign({:ep => :ep, :em => :em})

# e+e- selection: P > 0.7*E_beam, E > 0.6*P
# (applied in ROOT via algorithm.note)
alg_ee.with_decay_card(decay_card_ee).apply(event_selection_ee)

alg_ee.note(:ee_selection,
  "e+e- -> e+e- selection: |cos theta|<0.8, Vz<10 cm, Vr<1 cm. Both tracks: P > 0.7*E_beam, E > 0.6*P. Efficiency ~70%, background <0.05%. Cross sections measured at 15 energy points. Applied in ROOT.")

alg_ee.note(:cross_section_ee,
  "Cross section measurement: sigma = (Nsig-Nbkgs)/(L*eps_trg*eps_recsel)*f. Trigger efficiency ~100%. Efficiency correction factor f ~0.999. Measured cross sections 211-368 nb across energy points. Applied in ROOT.")

### Event selection: e+e- -> mu+mu- channel ###
alg_mumu = Algorithm.new("Jpsi_scan_mumu")
alg_mumu.set_header(["JpsiScanMuMuAlg/JpsiScanMuMu.h"])

event_selection_mumu = Selection.new

# Two oppositely charged tracks, barrel region
event_selection_mumu.select_track {
                      cos_theta 0.8
                      Vz   10.0
                      Vr   1.0
                      nTot "==2"
                    }
                    .assign({:mup => :mup, :mum => :mum})

# mu+mu- selection: P > 0.8*E_beam, 25 MeV < E < 0.25*P, no neutral showers > 25 MeV
# TOF timing: Delta_t < 1.5 ns (cosmic suppression)
# (applied in ROOT via algorithm.note)
alg_mumu.with_decay_card(decay_card_mumu).apply(event_selection_mumu)

alg_mumu.note(:mumu_selection,
  "mu+mu- selection: |cos theta|<0.8, Vz<10 cm, Vr<1 cm. Both tracks: P > 0.8*E_beam, 25 MeV < E < 0.25*P. No neutral showers E>25 MeV. |Delta_t_TOF| < 1.5 ns. Efficiency ~80%, background <0.5%. Cross sections measured at 15 energy points. Applied in ROOT.")

alg_mumu.note(:cross_section_mumu,
  "Cross section measurement: sigma = (Nsig-Nbkgs)/(L*eps_trg*eps_recsel)*f. Trigger efficiency ~100%. Efficiency correction factor f ~1.009. Measured cross sections 8-218 nb across energy points. Applied in ROOT.")

# Global chi2 fit: simultaneous fit to ee and mumu cross sections
alg_ee.note(:global_fit,
  "Global chi2 fit: simultaneous fit to cross sections of e+e- -> e+e- and e+e- -> mu+mu- at 15 energy points. Covariance matrices include correlations. Six free parameters: M, Gamma_ee*Gamma_ee/Gamma_tot, Gamma_ee*Gamma_mumu/Gamma_tot, S_W, F_ee, F_mumu. chi2_min/ndf = 23.0/24. Results: Gamma_ee*Gamma_ee/Gamma_tot = 0.346+/-0.009 keV, Gamma_ee*Gamma_mumu/Gamma_tot = 0.335+/-0.006 keV. Applied in ROOT.")

# 15 energy points
alg_ee.note(:energy_points,
  "15 CM energy points: 3049.642, 3058.693, 3082.496, 3088.854, 3091.760, 3094.697, 3095.430, 3095.826, 3097.213, 3098.340, 3099.042, 3101.359, 3105.580, 3112.051, 3119.878 MeV. Total integrated luminosity ~83 pb^-1. CM energies measured by BEMS, calibrated to J/psi PDG mass. Applied in ROOT.")

all_datasets_ee = all_incMC + exMC_ee
all_datasets_mumu = all_incMC + exMC_mumu
root_files_ee = alg_ee.execute_on(all_datasets_ee)
root_files_mumu = alg_mumu.execute_on(all_datasets_mumu)