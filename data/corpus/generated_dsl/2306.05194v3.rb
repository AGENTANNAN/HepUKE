# frozen_string_literal: true

### Dataset description ###
# Real data: eight c.m. energy points from 4.128 to 4.226 GeV (total ~7.33 fb^-1)
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

# Inclusive MC at each energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card: Ds+ -> eta e+ nu_e (eta -> gamma gamma); the opposite Ds- tag side uses hadronic modes
decay_card_eta = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 eta e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: Ds+ -> eta' e+ nu_e (eta' -> gamma rho0, rho0 -> pi+ pi-)
decay_card_etap = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s+ D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 eta' e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay eta'
  1.000 gamma rho0 PHSP;
  Enddecay

  Decay rho0
  1.000 pi+ pi- VSS;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive MC for each signal channel, one sample per energy point
exMC_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsToEtaENu"
  config.events        = 500000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

exMC_etap = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "DsToEtaPrimeENu"
  config.events        = 500000
  config.decay_card    = decay_card_etap
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg = TagAnalysis.new("DsToEtaENu")
alg.set_header(["DsToEtaENuAlg/DsToEtaENu.h"])
   .set_constant({"ECMS" => [:double, 4.178]})
   .with_decay_card(decay_card_eta)

# Ds- tag side: 14 hadronic tag modes (single tag; the tag carries its own
# selected, PID'd tracks and showers)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKPiPi, :DstoPiPiPi, :DstoKKPiPi0,
          :DstoKPiPiPi0, :DstoPiPiPiPi0, :DstoKsK, :DstoKsKPi0,
          :DstoKsPiPiPi, :DstoKsPiPiPiPi0, :DstoEtaPi, :DstoEtaPrimePi,
          :DstoPhiPi, :DstoPhiPiPi0
  t.charm -1   # pin the tagged anti-Ds (Ds-)
end

# Signal side: one positron, photons from eta -> gamma gamma, and a missing neutrino
alg.signal_side do |s|
  s.charged(ep: 1)          # exactly one positron
  s.photons 2               # eta -> gamma gamma (fit consumes up to 2 photons)
  s.min_photon_angle 10.0   # photon isolation from charged tracks
  s.missing :nu_e           # missing neutrino (massless)
end

# 3C kinematic fit: energy-momentum conservation + Ds mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures captured for the systematics stage
alg.note(:pid_correction_method,
  "signal-side positron identification requires CL_e > 0.001 and " \
  "CL_e/(CL_e+CL_pi+CL_K) > 0.8; the TagAnalysis signal-side lepton PID uses " \
  "fixed SimplePIDSvc thresholds in DSL v1")
alg.note(:background_veto,
  "extra-energy and extra-track vetoes applied on the signal side to reject " \
  "mis-reconstructed and continuum backgrounds")
alg.note(:efficiency_curve,
  "the transition gamma/pi0 is selected as the candidate giving the smallest |deltaE|; " \
  "the eta' -> gamma rho0 mode additionally contributes a third (transition) photon")
alg.note(:tag_track_selection,
  "tag-side charged tracks required to have |cos(theta)|<0.93, |Vz|<10 cm, " \
  "Vxy<1 cm; kaon/pion separation by likelihood ratios L_K>L_pi (kaons) and " \
  "L_pi>L_K (pions), applied inside DTagAlg")
alg.note(:tag_intermediate_reconstruction,
  "K_S0 -> pi+pi- accepted within |M-M(K_S0)|<12 MeV/c2; pi0/eta -> gamma gamma " \
  "mass windows applied inside DTagAlg")
alg.note(:mBC_window,
  "single-tag yields are extracted from M_BC cuts per energy point and M_tag fits; " \
  "the tag mBC is stored unconditionally and windowed in the ROOT fit")
alg.note(:post_fit_cut,
  "M(eta(')e+) < 1.9 GeV/c2 is applied after the 3C kinematic fit and evaluated in ROOT")
alg.note(:signal_extraction,
  "the signal is extracted with a simultaneous M_miss^2 fit over the eta and eta' channels in ROOT")

alg.apply
alg.execute_on(data_points + incMC_points + exMC_eta + exMC_etap)