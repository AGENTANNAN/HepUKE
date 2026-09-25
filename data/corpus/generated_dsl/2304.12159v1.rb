# Core DSL classes and dependencies are loaded automatically at execution time.

### Dataset description ###
# Real data at eight c.m. energies (4.128 - 4.226 GeV) spanning the psi(4260) region
data_4130 = DatasetManager.real_data.find("705_4130")   # 4.128 GeV
data_4160 = DatasetManager.real_data.find("705_4160")   # 4.157 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV (psi(4260) peak)

# Corresponding inclusive MC samples
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

# Decay card for the signal process (EvtGen format):
#   psi(4260) -> D_s*+ D_s- , D_s*+ -> e+ nu_e , D_s- -> K+ K- pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0 e+ nu_e PHSP;
    Enddecay

    Decay D_s-
    1.0 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Single 500k-event exclusive signal MC sample, related to the psi(4260) peak point
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4260_dsstar_enu"
  config.related_dataset = data_4230
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (TagAnalysis) ###
alg_name = "DsStarENuTag"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.226]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .with_decay_card(decay_card_signal)

# Tag side: hadronic D_s- tag (charm = -1).
# A single tag_side declaration yields the single-tag (ST) pattern.
# Nine of the nominal 16-mode ST list are applied; the undeclared modes are
# flagged unavailable (kept out of the reconstruction mode lists).
tag_alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,                 # K+ K- pi-
          :DstoKKPiPi0,              # K+ K- pi- pi0
          :DstoPiPiPi,               # pi+ pi- pi-
          :DstoKsK,                  # K_S0 K-
          :DstoKsKPiPi,              # K_S0 K+ pi- pi-
          :DstoPiEta,                # pi- eta
          :DstoRhoEta,               # rho- eta
          :DstoPiEtaPrime,           # pi- eta'(-> pi+ pi- eta)
          :DstoPiEtaPrimeGammaRho    # pi- eta'(-> gamma rho0)
  t.charm -1
end

# Signal side: exactly one extra charged track identified as e+, charge +1,
# and one massless missing nu_e (semileptonic signal opposite the D_s- tag).
tag_alg.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
  s.min_photon_angle 10.0
end

# 4C kinematic fit against the measured CMS four-momentum, chi2 < 200
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side / DTagAlg-internal procedures with no corresponding DSL method
tag_alg
  .note(:tag_track_selection, "ST charged tracks: |cos(theta)| < 0.93, |Vxy| < 1 cm, |Vz| < 10 cm (20 cm for K_S0 daughters); applied inside DTagAlg")
  .note(:pid_correction_method, "K/pi separation in the D_s- tag uses the likelihood ratios L(K/pi) > L(pi/K) and L(K/pi) > L(e); implemented inside DTagAlg")
  .note(:secondary_particle_reconstruction, "K_S0 from pi+pi- with |M - M_K_S0| < 12 MeV/c2 and decay length > 2 sigma; pi0/eta from gamma-gamma with mass windows 115-150 (pi0) / 500-570 (eta) MeV/c2 plus a mass-constrained fit; rho mass 570-970 MeV/c2; eta' windows M(eta pi+pi-) = 946-970 and M(gamma rho0) = 940-976 MeV/c2; all reconstructed inside DTagAlg")
  .note(:background_veto, "K_S0 veto applied in the pi+pi-pi- and K-pi+pi- tag modes")
  .note(:tag_candidate_selection, "direct pions from D_s- required to have p > 100 MeV/c; the tag candidate closest to the PDG D_s- mass is kept, with M(D_s-) within 3 sigma")
  .note(:signal_electron_pid, "DT signal electron: exactly one extra charged track with L_e > 0.8*(L_e + L_pi + L_K) and L_e > 0.001; bremsstrahlung recovery adds EMC showers within 10 deg of the e+; lepton PID thresholds are fixed in the generated DTag code")
  .note(:efficiency_curve, "M_rec signal region to be defined by S/sqrt(S+B) optimisation (applied at the ROOT analysis stage)")

# Render the tag specification (zero-argument apply — no Selection object)
tag_alg.apply

# Execute on all eight real-data points, their inclusive MC, and the signal MC
root_files = tag_alg.execute_on([
  data_4130, data_4160, data_4180, data_4190,
  data_4200, data_4210, data_4220, data_4230,
  incMC_4130, incMC_4160, incMC_4180, incMC_4190,
  incMC_4200, incMC_4210, incMC_4220, incMC_4230,
  exMC_signal
])