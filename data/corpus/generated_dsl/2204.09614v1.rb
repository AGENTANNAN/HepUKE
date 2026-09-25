# =============================================================================
# e+e- -> D_s*+ D_s- -> gamma D_s+ D_s-,  D_s+ -> K_S0 K+ pi0
# (K_S0 -> pi+ pi-, pi0 -> gamma gamma)  at sqrt(s) = 4.178 - 4.226 GeV
# Tag-based BOSS spec: the D_s- is tagged through seven hadronic modes, the
# signal D_s+ is reconstructed on the side the tag did not use, and the D_s*+
# transition photon accompanies the signal.
# =============================================================================

### Dataset preparation ###
data_names = ["703_4180", "703_4190", "703_4200", "703_4210", "703_4220", "703_4230"]
data  = data_names.map { |n| DatasetManager.real_data.find(n) }     # six energy points
incmc = data_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # matching inclusive MC

# Decay card for the full signal chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.000 K_S0 K+ pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Signal MC for every energy point (same decay chain, 2M events in total)
exMCs = DatasetManager.create_exclusive_mc_for(data) do |config|
  config.sample_name   = "exmc_Dsstar_Ds_gamma_KsKpi"
  config.events        = 2_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag analysis (BOSS) ###
alg = TagAnalysis.new("DsTagKsKpi")
alg.set_header(["DsTagKsKpiAlg/DsTagKsKpi.h"])
   .set_constant({ "ECMS" => [:double, 4.226] })   # measured beam energy is used on data; ECMS is the MC fallback

# ---- Tag side: D_s- from the pre-stored DTagAlg candidates ------------------
# Seven hadronic modes; the eighth mode (K_S0 K+ pi- pi-) is not available.
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoPiPiPi, :DstoEtaPrimePi, :DstoKPiPi
  t.charm(-1)                              # pin the tagged side to D_s-
  t.window :mBC, min: 1.930, max: 1.990    # explicit tag-candidate M_tag window
end

# ---- Signal side: tracks / showers not used by the tag ----------------------
alg.signal_side do |s|
  s.photons 2                              # two pi0 photons
  s.charged(kp: 1, pip: 1, pim: 1)         # K+, pi+ (from K_S0), pi-
  s.require_charge(1)                      # net signal-side charge = +1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# ---- Kinematic fit: 4C + mass constraints -----------------------------------
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)                 # K_S0 constraint
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)              # pi0 constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)                        # D_s- tag constraint
  f.invariant_mass_of(:pip, :pim, :kp, :gamma, :gamma).constrain_to_nominal_mass_of(:Ds)  # added D_s+ constraint
  f.chi2_cut 200
end

# ---- BOSS-side procedures that the tag DSL cannot express -------------------
alg.note(:ks0_selection,
         "K_S0 from pi+ pi-: daughter pion |Vz| < 20 cm and |cos(theta)| < 0.93; " \
         "secondary-vertex fit chi2 < 100; K_S0 flight significance > 2 sigma; " \
         "pi+ pi- invariant mass in (0.492, 0.503) GeV/c^2")
alg.note(:pi0_selection,
         "pi0 from gamma gamma through a 1C gamma-gamma mass-constrained fit with chi2 < 10; " \
         "photons with E > 25 MeV (barrel, |cos theta| < 0.80) or E > 50 MeV " \
         "(endcap, 0.86 < |cos theta| < 0.92) and EMC timing < 700 ns")
alg.note(:dsstar_constraint,
         "the D_s*+ transition photon is not a fit participant (v1 consumes at most two " \
         "signal-side photons and both are taken by the pi0), so the D_s*+ mass constraint - " \
         "the 8th constraint of the 8C fit - is applied in ROOT against the recoil of the D_s- tag")
alg.note(:background_veto,
         "combinatoric suppression: reject M(K+ K-) < 1.05 GeV (tag side); " \
         "|M(K+ pi-) - M(K*(892))| < 70 MeV; |M(pi- pi0) - M(rho)| < 150 MeV")
alg.note(:tag_candidate_selection,
         "single tag: choose the candidate whose recoil mass is closest to the D_s*+ mass; " \
         "double tag: choose the candidate whose (M_tag + M_sig)/2 is closest to the D_s mass; " \
         "tag mBC/deltaE are stored unconditionally and this selection is performed in ROOT")
alg.note(:efficiency_curve,
         "D_s+ momentum > 0.1 GeV/c required to reject soft pions from D*+ decays")

alg.apply
alg.execute_on(data + incmc + exMCs)