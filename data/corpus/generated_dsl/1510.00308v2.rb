# ================================================================
# e+e- -> D+ D- at psi(3770) (3.773 GeV)
#   tag side : hadronic single-tag D (both charges scanned)
#   signal   : semileptonic D -> K_L0 e+ nu_e (K_L0 left missing)
# ================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 2.92 fb^-1 psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

### Decay cards + exclusive MC (500k events per tag mode) ###
# Tag-side hadronic modes; DTagTool scans both tag charges.
tag_mode_decays = {
  "DptoKPiPi"    => "K+ pi- pi-",       # K+ pi- pi-
  "DptoKPiPiPi0" => "K+ pi- pi- pi0",   # K+ pi- pi- pi0
  "DptoKsPiPi0"  => "K_S0 pi- pi0",     # K_S0 pi- pi0
  "DptoKsPiPiPi" => "K_S0 pi- pi- pi+", # K_S0 pi- pi- pi+
  "DptoKsPi"     => "K_S0 pi-",         # K_S0 pi-
  "DptoKKPi"     => "K+ K- pi-"         # K+ K- pi-
}

exMC_signals = tag_mode_decays.map do |mode_name, tag_decay|
  # e+e- -> D+ D-, D+ -> K_L0 e+ nu_e (signal), D- -> hadronic tag mode
  decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 K_L0 e+ nu_e PHSP;
    Enddecay

    Decay D-
    1.0 #{tag_decay} PHSP;
    Enddecay

    End
  DECAYCARD

  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_#{mode_name}_Kslnu"
    config.related_dataset = data_3773
    config.events          = 500_000
    config.decay_card      = decay_card
    config.cross_section   = :default
  end
end

### Tag-based event selection (BOSS) ###
alg_name = "DpTagKslnu"
d_tag_alg = TagAnalysis.new(alg_name)
d_tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })

# Hadronic single-tag D: the six modes, both charges scanned (no charm pin)
d_tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKsPi, :DptoKKPi
  # Explicit tag-side windows requested by the analysis (default would be store-not-cut)
  t.window :mBC,    min: 1.86, max: 1.88          # 1.86 < M_BC < 1.88 GeV
  t.window :deltaE, min: -0.057, max: 0.040       # envelope of the mode-dependent ~+-3 sigma windows
end

# Signal side: D+ -> K_L0 e+ nu_e
d_tag_alg.signal_side do |s|
  s.charged(ep: 1)         # exactly one positron, charge opposite the tag D; no other unused track
  s.photons 1..1           # exactly one unused photon (defines the K_L0 direction)
  s.min_photon_energy 0.1  # E_gamma > 0.1 GeV
  s.min_photon_angle 10.0  # photon angle to charged tracks > 10 degrees
  s.missing :K_L0          # K_L0 left missing; |p| solved from 4-momentum conservation (U_miss = 0)
end

# Kinematic fit: four-momentum constraint, chi2 < 200
d_tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures
d_tag_alg
  .note(:mode_dependent_deltaE,
        "Tag selection uses mode-dependent deltaE windows (~+-3 sigma): K pi pi +-30 MeV, " \
        "K pi pi pi0 -52/+39 MeV, Ks pi pi0 -57/+40 MeV, Ks pi pi pi +-34 MeV, Ks pi +-32 MeV, " \
        "K K pi +-30 MeV; per mode only the candidate with the smallest |deltaE| is kept.")
  .note(:tag_pi0_Ks_reconstruction,
        "Tag pi0: photon pairs with 0.110 < M(gamma gamma) < 0.155 GeV plus a 1C " \
        "mass-constrained fit to the nominal pi0 mass with chi2 < 20. Tag K_S0: oppositely " \
        "charged track pairs with a vertex fit, vertex within 20 cm along the beam axis, " \
        "|M(pi+pi-) - m(K_S0)| < 12 MeV and decay length greater than 2 sigma from the IP.")
  .note(:electron_id,
        "Signal positron identified by the combined EMC/dE/dx/TOF likelihood ratio " \
        "L'(e)/[L'(e)+L'(pi)+L'(K)] > 0.8, with bremsstrahlung recovery adding unmatched " \
        "showers within 5 degrees of the electron.")
  .note(:background_veto,
        "Signal photons pairing into 0.110 < M(gamma gamma) < 0.155 GeV are vetoed to suppress " \
        "pi0 feed-through; the most energetic qualifying shower defines the K_L0 direction.")

d_tag_alg.apply
root_files = d_tag_alg.execute_on([data_3773, incMC_3773] + exMC_signals)