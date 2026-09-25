# ============================================================================
# Datasets: eight c.m. energy points, 4.128 - 4.226 GeV
# ============================================================================
data_points = [
  DatasetManager.real_data.find("705_4130"),   # 4.128 GeV
  DatasetManager.real_data.find("705_4160"),   # 4.157 GeV
  DatasetManager.real_data.find("703_4180"),   # 4.178 GeV
  DatasetManager.real_data.find("703_4190"),   # 4.189 GeV
  DatasetManager.real_data.find("703_4200"),   # 4.199 GeV
  DatasetManager.real_data.find("703_4210"),   # 4.209 GeV
  DatasetManager.real_data.find("703_4220"),   # 4.219 GeV
  DatasetManager.real_data.find("703_4230")    # 4.226 GeV
]

inc_mc_points = [
  DatasetManager.inclusive_mc.find("705_4130"),  # 4.128 GeV
  DatasetManager.inclusive_mc.find("705_4160"),  # 4.157 GeV
  DatasetManager.inclusive_mc.find("703_4180"),  # 4.178 GeV
  DatasetManager.inclusive_mc.find("703_4190"),  # 4.189 GeV
  DatasetManager.inclusive_mc.find("703_4200"),  # 4.199 GeV
  DatasetManager.inclusive_mc.find("703_4210"),  # 4.209 GeV
  DatasetManager.inclusive_mc.find("703_4220"),  # 4.219 GeV
  DatasetManager.inclusive_mc.find("703_4230")   # 4.226 GeV
]

# ============================================================================
# Signal decay cards (EvtGen format)
# ============================================================================
# Mode I: e+e- -> Ds*+ Ds- ; Ds*+ -> gamma Ds+ ; Ds+ -> f1(1285) e+ nu_e ;
#         f1(1285) -> pi+ pi- eta ; eta -> gamma gamma
decay_card_f1285 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.0000 f_1(1285) e+ nu_e PHSP;
  Enddecay

  Decay f_1(1285)
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: e+e- -> Ds*+ Ds- ; Ds*+ -> gamma Ds+ ; Ds+ -> f1(1420) e+ nu_e ;
#          f1(1420) -> K+ K- pi0 ; pi0 -> gamma gamma
decay_card_f1420 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.0000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.0000 f_1(1420) e+ nu_e PHSP;
  Enddecay

  Decay f_1(1420)
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# ============================================================================
# Exclusive MC: 500k events per signal mode, one sample per energy point
# ============================================================================
exMC_f1285 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ds_f1285_ep_nu_exmc"
  config.events        = 500_000
  config.decay_card    = decay_card_f1285
  config.cross_section = :default
end

exMC_f1420 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ds_f1420_ep_nu_exmc"
  config.events        = 500_000
  config.decay_card    = decay_card_f1420
  config.cross_section = :default
end

# ============================================================================
# Chain I : Ds+ -> f1(1285) e+ nu_e, f1(1285) -> pi+ pi- eta, eta -> gamma gamma
# ============================================================================
alg_f1285 = TagAnalysis.new("DsTagF1285")
alg_f1285.set_header(["DsTagF1285Alg/DsTagF1285.h"])
         .set_constant({"ECMS" => [:double, 4.178]})
         # >= 3 signal-side photons are required (2 from eta -> gamma gamma plus the
         # transition photon of Ds*+- -> gamma Ds+-); the tag-fit participant list
         # consumes at most two photons, so the transition photon is kept at the
         # event level and used in the ROOT stage.
         .note(:tag_transition_photon, "signal side requires a third photon, the transition photon of Ds*+- -> gamma Ds+- (Egamma > 25 MeV, angle to the nearest charged track > 10 deg); the tag kinematic fit consumes at most two photons, so the transition photon and its Ds*+- parent assignment are resolved in the ROOT stage")
         # Electron identification criteria (combined MDC+TOF+EMC probabilities,
         # L'_e > 0, L'_e/(L'_e+L'_pi+L'_K) > 0.8, E_EMC/|p| > 0.8) are not tunable
         # through the tag signal-side electron recipe (fixed v1 thresholds).
         .note(:electron_pid, "signal e+ identified by combined MDC+TOF+EMC probabilities with L'_e > 0, L'_e/(L'_e+L'_pi+L'_K) > 0.8 and E_EMC/|p| > 0.8 c; the tag signal-side electron selection uses the fixed v1 SimplePIDSvc recipe instead of these dedicated cuts")
         # Non-K_S0 pion track quality.
         .note(:non_ks_pion_selection, "signal-side pi+ and pi- (non-K_S0) required to have momentum p > 100 MeV/c")
         # Per-energy missing-mass windows: stored by the tag fit and applied in ROOT.
         .note(:missing_mass_window, "missing mass against the tagged Ds- (M_rec) is stored and windowed per energy point in the ROOT stage, e.g. [2.048, 2.190] GeV/c^2 at 4.178 GeV; the M_rec^2 variable built from the tagged Ds- plus the transition photon is windowed at (3.82, 4.05) GeV^2/c^4")
         # Two Ds* parent-assignment hypotheses.
         .note(:ds_star_parent_hypothesis, "the transition photon is assigned either to Ds*+ -> gamma Ds+ (signal side) or to Ds*- -> gamma Ds- (tag side); the hypothesis giving the smallest 3C chi2 of the fit that also constrains the Ds-, Ds+ and Ds* masses is retained")
         .with_decay_card(decay_card_f1285)

# Ds- tag: twelve hadronic tag modes (hadronic tag built from the pre-stored tag
# candidates of the DTagTool; the Ds- side is pinned by charm = -1)
alg_f1285.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoPiEta, :DstoPiEtaPrime, :DstoKKPiPi0,
          :DstoKsKPiPi, :DstoRhoEta, :DstoPiEtaPrimeRho, :DstoKsKPi0,
          :DstoPiPiPi, :DstoKPiPi
  t.charm -1
end

# Signal side: everything the tag did not use -> e+ pi+ pi- + eta(gamma gamma) + nu_e
alg_f1285.signal_side do |s|
  s.charged(ep: 1, pip: 1, pim: 1)   # one e+, one pi+, one pi- (f1(1285) channel)
  s.photons 2                        # two photons forming the eta -> gamma gamma
  s.min_photon_angle 10.0            # photon angle to the nearest charged track > 10 deg
  s.min_photon_energy 0.025          # Egamma > 25 MeV
  s.missing :nu_e                    # semileptonic decay -> one missing electron neutrino
  s.require_charge 1                 # net charge of the signal side = +1
end

# Kinematic fit: 3C (four-momentum conservation with a missing nu_e) plus the
# eta -> gamma gamma mass constraint; also constraining the tagged Ds- mass and
# the signal Ds+ mass.
alg_f1285.fit do |f|
  f.constrain_four_momentum                                                   # 4-momentum conservation (neutrino floated)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)      # eta -> gamma gamma
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)                # tagged Ds- mass
  f.invariant_mass_of(:ep, :pip, :pim, :gamma, :gamma, :nu_e)
   .constrain_to_nominal_mass_of(:Ds)                                         # signal Ds+ mass
  f.chi2_cut 200
  f.store_fitted_momenta                                                      # keep fitted four-momenta
end

alg_f1285.apply
root_files_f1285 = alg_f1285.execute_on(data_points + inc_mc_points + exMC_f1285)

# ============================================================================
# Chain II : Ds+ -> f1(1420) e+ nu_e, f1(1420) -> K+ K- pi0, pi0 -> gamma gamma
# ============================================================================
alg_f1420 = TagAnalysis.new("DsTagF1420")
alg_f1420.set_header(["DsTagF1420Alg/DsTagF1420.h"])
         .set_constant({"ECMS" => [:double, 4.178]})
         .note(:tag_transition_photon, "signal side requires a third photon, the transition photon of Ds*+- -> gamma Ds+- (Egamma > 25 MeV, angle to the nearest charged track > 10 deg); the tag kinematic fit consumes at most two photons, so the transition photon and its Ds*+- parent assignment are resolved in the ROOT stage")
         .note(:electron_pid, "signal e+ identified by combined MDC+TOF+EMC probabilities with L'_e > 0, L'_e/(L'_e+L'_pi+L'_K) > 0.8 and E_EMC/|p| > 0.8 c; the tag signal-side electron selection uses the fixed v1 SimplePIDSvc recipe instead of these dedicated cuts")
         .note(:non_ks_pion_selection, "signal-side K+ and K- required to have momentum p > 100 MeV/c")
         .note(:missing_mass_window, "missing mass against the tagged Ds- (M_rec) is stored and windowed per energy point in the ROOT stage; the M_rec^2 variable built from the tagged Ds- plus the transition photon is windowed at (3.78, 4.05) GeV^2/c^4")
         .note(:ds_star_parent_hypothesis, "the transition photon is assigned either to Ds*+ -> gamma Ds+ (signal side) or to Ds*- -> gamma Ds- (tag side); the hypothesis giving the smallest 3C chi2 of the fit that also constrains the Ds-, Ds+ and Ds* masses is retained")
         .note(:background_veto, "Ds+ -> phi e+ nu_e suppressed by requiring M(K+K-) > 1.03 GeV/c^2; Ds+ -> K+K-pi+pi0 background suppressed by requiring M(K+K-pi0 e+) > 1.92 GeV/c^2")
         .with_decay_card(decay_card_f1420)

# Ds- tag: nine of the hadronic tag modes (the two eta/eta' pion modes dropped)
alg_f1420.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKKPiPi0, :DstoKsKPiPi, :DstoRhoEta,
          :DstoPiEtaPrimeRho, :DstoKsKPi0, :DstoPiPiPi, :DstoKPiPi
  t.charm -1
end

# Signal side: e+ K+ K- + pi0(gamma gamma) + nu_e
alg_f1420.signal_side do |s|
  s.charged(ep: 1, kp: 1, km: 1)     # one e+, one K+, one K- (f1(1420) channel)
  s.photons 2                        # two photons forming the pi0 -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025          # Egamma > 25 MeV
  s.missing :nu_e                    # semileptonic decay -> one missing electron neutrino
  s.require_charge 1                 # net charge of the signal side = +1
end

# Kinematic fit: 4-momentum conservation + pi0 -> gamma gamma mass constraint,
# plus the tagged Ds- and the signal Ds+ mass constraints.
alg_f1420.fit do |f|
  f.constrain_four_momentum                                                   # 4-momentum conservation (neutrino floated)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)      # pi0 -> gamma gamma
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)                # tagged Ds- mass
  f.invariant_mass_of(:ep, :kp, :km, :gamma, :gamma, :nu_e)
   .constrain_to_nominal_mass_of(:Ds)                                         # signal Ds+ mass
  f.chi2_cut 200
  f.store_fitted_momenta                                                      # keep fitted four-momenta
end

alg_f1420.apply
root_files_f1420 = alg_f1420.execute_on(data_points + inc_mc_points + exMC_f1420)