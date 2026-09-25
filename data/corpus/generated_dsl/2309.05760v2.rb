# ============================================================================
# ψ(3770) at 3.773 GeV : e+e- -> D+ D-
# Signal      : D+ -> K_S0 π+ η,  K_S0 -> π+π-,  η -> γγ
# Double tag  : the recoiling D- is reconstructed from the pre-stored DTag
#               candidates in six hadronic modes (the tag supplies its own
#               selected, PID'ed tracks and showers)
# ============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data (BOSS 712)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC sample

# Decay card for the exclusive signal MC of the tagged process
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0000 K_S0 pi+ eta PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC for D+ -> K_S0 π+ η
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DplusToKsPiEta"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DplusToKsPiEta"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770) CMS energy in GeV
   .with_decay_card(decay_card_signal)

# --- Tag side : the recoiling D- in six hadronic modes (charm -1 pins the D-) ---
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,      # D- -> K+ π- π-
          :DptoKsPi,       # D- -> K_S0 π-
          :DptoKPiPiPi0,   # D- -> K+ π- π- π0
          :DptoKsPiPi0,    # D- -> K_S0 π- π0
          :DptoKsPiPiPi,   # D- -> K_S0 π- π- π+
          :DptoKKPi        # D- -> K+ K- π-
  t.charm -1                                   # tag the D- (negative charm side)
  t.window :mBC, min: 1.865, max: 1.875        # M_BC ∈ [1.865, 1.875] GeV (uniform for all tag modes)
end

# --- Signal side : D+ -> K_S0 π+ η with K_S0 -> π+π-, η -> γγ ---
alg.signal_side do |s|
  s.photons 2                 # the two photons from η -> γγ
  s.min_photon_angle 10.0     # photon opening angle > 10 degrees
  s.min_photon_energy 0.025   # photon energy > 25 MeV
  s.charged(pip: 2, pim: 1)   # π+ (from D+) + π+π- (from K_S0)
  s.require_charge 1          # net charge +1
end

# --- 4C kinematic fit : γγ → η, π+π- → K_S0, tag D- → nominal D mass ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D-")
  f.chi2_cut 200
end

# --- BOSS-side procedures that cannot be expressed in the tag DSL ---
alg
  .note(:tag_deltaE_window, "tag ΔE is stored unconditionally; the mode-dependent window " \
        "([-0.025, 0.025] GeV for the non-π0 tag modes K+π-π-, K_S0π-, K_S0π-π-π+, K+K-π-; " \
        "[-0.055, 0.040] GeV for the π0 tag modes K+π-π-π0, K_S0π-π0) is applied in ROOT, " \
        "because only one window per observable is allowed on a single tag side")
  .note(:signal_deltaE_window, "signal-side ΔE window [-0.020, 0.020] GeV applied in ROOT on " \
        "the stored signal-side ΔE")
  .note(:best_candidate_selection, "best tag and signal candidates chosen by minimum |ΔE|; " \
        "DTag ranking supports only :inv / :mbc and the fit selects the signal combination by " \
        "χ², so the minimum-|ΔE| choice is performed in ROOT")
  .note(:background_veto, "events containing an extra EMC shower with energy above 0.23 GeV " \
        "(not among the two signal photons) are vetoed")

# Validate + render the tag specification (no Selection argument), then run
alg.apply
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])