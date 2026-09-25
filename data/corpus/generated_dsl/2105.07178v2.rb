### Dataset preparation ###
# Real data at the six energy points of the e+e- -> Ds+ Ds- scan
data_4178 = DatasetManager.real_data.find("703_4180")   # 4178 MeV  (4.178 GeV)
data_4189 = DatasetManager.real_data.find("703_4190")   # 4188.8 MeV (4.189 GeV)
data_4199 = DatasetManager.real_data.find("703_4200")   # 4198.9 MeV (4.199 GeV)
data_4209 = DatasetManager.real_data.find("703_4210")   # 4209.2 MeV (4.209 GeV)
data_4219 = DatasetManager.real_data.find("703_4220")   # 4218.7 MeV (4.219 GeV)
data_4226 = DatasetManager.real_data.find("703_4230")   # 4226.3 MeV (4.226 GeV)

# Inclusive MC at 4.178 GeV
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card for the signal chain:
#   e+e- -> Ds+ Ds-,  Ds+ -> tau+ nu_tau,  tau+ -> pi+ pi0 anti-nu_tau,  pi0 -> gamma gamma
# (top mother psi(4260) per BESIII/KKMC convention; the tag side Ds- is generated
#  through a representative hadronic channel, the full tag-mode list is used at
#  reconstruction time via DTagTool).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Ds+ Ds- PHSP;
    Enddecay

    Decay Ds+
    1.0 tau+ nu_tau PHSP;
    Enddecay

    Decay tau+
    1.0 pi+ pi0 anti-nu_tau PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay Ds-
    1.0 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 2,000,000 exclusive MC events for the full anti-Ds- hadronic-tag + Ds+ -> tau+ nu_tau chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_ds_taunu_pipipi0nu"
  config.related_dataset = data_4178      # anchored at the 4.178 GeV reference point
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (tag-based, BOSS) ###
alg_name = "DsTaunuTag"
ds_taunu = TagAnalysis.new(alg_name)
ds_taunu.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 4.178]})   # reference energy point for the fit setup
        .with_decay_card(decay_card_signal)

# Tag side: hadronic Ds- reconstructed from pre-stored tag candidates (single tag),
# 14 DTagAlg hadronic modes.
ds_taunu.tag_side(:Ds) do |t|
  t.modes :DstoKsK,          # K_S0 K-
          :DstoKsKPi0,       # K_S0 K- pi0
          :DstoKsKPiPi,      # K_S0 K- pi+ pi-  /  K_S0 K+ pi- pi-
          :DstoKsKPiPiPi,    # K_S0 K- pi+ pi- pi0
          :DstoKKPi,         # K+ K- pi-
          :DstoKKPiPi0,      # K+ K- pi- pi0
          :DstoKKPiPiPi,     # K+ K- pi- pi+ pi-
          :DstoPiPiPi,       # pi+ pi- pi-
          :DstoPiPiPiPi0,    # pi+ pi- pi- pi0
          :DstoPiPiPiPiPi,   # pi+ pi+ pi- pi- pi-
          :DstoKsPi,         # K_S0 pi-
          :DstoKsPiPi0,      # K_S0 pi- pi0
          :DstoKsPiPiPi,     # K_S0 pi- pi+ pi-
          :DstoKsPiPiPiPi    # K_S0 pi- pi+ pi- pi0
  t.charm -1                 # pin the tagged side to Ds- (charm = -1)
end

# Signal side: everything the tag did not use -> Ds+ -> tau+ nu_tau, tau+ -> pi+ pi0 nu_tau
ds_taunu.signal_side do |s|
  s.charged(pip: 1)         # exactly one positive charged pion
  s.photons 2               # exactly two photons (from pi0 -> gamma gamma)
  s.min_photon_angle 10.0   # photon-track isolation angle (degrees)
  s.require_charge 1        # net charge +1
  s.missing :nu_tau         # one missing massless tau neutrino
end

# Kinematic fit: 4-momentum conservation + pi0 mass constraint on the two photons.
# MM^2 discrimination, extra-photon vetoes, pi0 mass window and the final chi2
# optimisation are deferred to the ROOT analysis (simultaneous fits over the six points).
ds_taunu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Render the tag analysis (takes no Selection argument)
ds_taunu.apply

# Execute on the six data points, the 4.178 GeV inclusive MC and the signal exclusive MC
root_files = ds_taunu.execute_on([data_4178, data_4189, data_4199,
                                  data_4209, data_4219, data_4226,
                                  incMC_4178, exMC_signal])