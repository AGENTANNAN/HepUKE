# ============================================================================
# Datasets: 20.3 fb^-1 psi(3770) data at sqrt(s) = 3.773 GeV + inclusive MC
# ============================================================================
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# ============================================================================
# Signal decay card (EvtGen):
#   psi(3770) -> D+ D- ,  D+ -> a0(980)+ eta , a0(980)+ -> pi+ eta , eta -> gamma gamma
# ============================================================================
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0000 a0(980)+ eta PHSP;
  Enddecay

  Decay a0(980)+
  1.0000 pi+ eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC for the full D+ -> pi+ eta eta chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToPiEtaEta"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ============================================================================
# Tag analysis: tag the recoil D- in six hadronic modes,
#               signal side D+ -> pi+ eta eta (both eta -> gamma gamma)
# ============================================================================
alg_name = "DpToPiEtaEtaTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# --- Tag side: recoil D- in six hadronic tag modes --------------------------
alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi,        # D- -> K+ pi- pi-
          :DptoKPiPiPi0,     # D- -> K+ pi- pi- pi0
          :DptoKsPi,         # D- -> K_S0 pi-
          :DptoKsPiPi0,      # D- -> K_S0 pi- pi0
          :DptoKsPiPiPi,     # D- -> K_S0 pi- pi- pi+
          :DptoKKPi          # D- -> K+ K- pi-
  t.charm -1                 # pin the tagged side to D-

  # Explicit pre-fit tag-side windows requested in the description
  t.window :mBC,    min: 1.83     # reject tag candidates with M_BC < 1.83 GeV/c^2
  t.window :deltaE, abs: 0.1      # reject tag candidates with |deltaE| > 0.1 GeV
end

# --- Signal side: exactly one pi+ and four photons --------------------------
alg.signal_side do |s|
  s.charged(pip: 1)          # exactly one positively charged pion
  s.photons 4                # four photons (2 x eta -> gamma gamma)
  s.min_photon_angle 10.0    # minimum photon opening angle of 10 degrees
end

# --- Kinematic fit: 4C + eta -> gamma gamma nominal-mass constraints ---------
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # eta -> gamma gamma (both candidates)
  f.chi2_cut 200
end

# --- BOSS-side procedures with no dedicated DSL construct -------------------
alg.note(:tag_eta_mass_window,
         "Tag-side eta reconstruction uses a widened mass window 0.45 < M(gamma gamma) < 0.65 GeV/c^2; " \
         "all other tag-side tracking, PID, K_S0 and pi0 handling follows the standard BESIII D-tag prescription.")
   .note(:tag_deltae_window_per_mode,
         "Per-mode +-3.5 sigma deltaE windows applied to the tag candidate before the kinematic fit; " \
         "the flat |deltaE| < 0.1 GeV and M_BC > 1.83 GeV/c^2 requirements are only the loose pre-selection.")
   .note(:best_tag_candidate,
         "One best tag candidate per tag mode is retained, chosen by |deltaE| closest to zero.")
   .note(:signal_prefit_window,
         "Signal-side pre-fit selection: 1.860 < M_BC < 1.880 GeV/c^2 and |deltaE| < 0.040 GeV; " \
         "the best signal candidate is kept by minimal |M_BC - m(D+)|, i.e. closest to the nominal D+ mass.")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])