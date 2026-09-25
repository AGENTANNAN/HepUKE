# ============================================================================
# ψ(3686) → π0 h_c analysis  (h_c → γ η′/η  and  h_c → γ π0 search)
# BOSS-side event selection up to (and including) the final kinematic fit.
# Five independent signal modes → five Algorithm objects (Rule T1).
# ============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # 27.12×10^8 ψ(3686) events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Corresponding inclusive MC

# ---------------------------------------------------------------------------
# Decay cards (EvtGen format) — one per signal mode
# ---------------------------------------------------------------------------
# Mode 1: h_c → γ η′, η′ → π+π−η, η → γγ
decay_card_mode1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.0000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 2: h_c → γ η′, η′ → γ π+π−
decay_card_mode2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.0000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 3: h_c → γ η, η → γγ
decay_card_mode3 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.0000 gamma eta PHSP;
    Enddecay
    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 4: h_c → γ η, η → π+π−π0
decay_card_mode4 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.0000 gamma eta PHSP;
    Enddecay
    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 5: h_c → γ π0 (search), π0 → γγ
decay_card_mode5 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.0000 gamma pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC: 100k events per mode
# ---------------------------------------------------------------------------
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_etapipiEta"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_etapigammapipi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_eta_gg"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_eta_pipipimpi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode4
  config.cross_section   = :default
end

exMC_mode5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_gammapi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode5
  config.cross_section   = :default
end

# ===========================================================================
### Event selection (BOSS) ###
# ===========================================================================

# ---------------------------------------------------------------------------
# Mode 1: η′ → π+π−η, η → γγ   (charged final state π+π−5γ)
#   final 6C fit (4C + π0 mass + η mass), χ² < 70
# ---------------------------------------------------------------------------
alg_name_mode1 = "HcEtaPipiEta"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode1 = Selection.new
  .select_track {                        # charged track selection
    cos_theta 0.93                       # |cosθ| < 0.93
    Vz        10.0                       # |Vz| < 10 cm
    Vr        1.0                        # Vr < 1 cm
    nChrp     "==1"                      # one positive track
    nChrn     "==1"                      # one negative track
    nNet      "==0"                      # net charge zero
  }
  .select_photon {                       # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0               # > 10° from nearest charged track
    energyThreshold_b 0.025              # 25 MeV (barrel)
    energyThreshold_e 0.050              # 50 MeV (endcap)
    nGam              ">=5"              # at least five photons
  }
  .pid(method: :probability) {
    prob_cut 0.001                       # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"                           # one π+
    npim "==1"                           # one π−
  }
  # intermediate Kalman fits: tag the η (γγ) and the π0 (γγ)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # final 6C kinematic fit (4C + π0 mass + η mass), loose χ² cut (tight cut in ROOT)
  .kinematic_fit([:pi0, :gamma, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 70
  }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.execute_on([psip_data, psip_incMC, exMC_mode1])

# ---------------------------------------------------------------------------
# Mode 2: η′ → γπ+π−   (charged final state π+π−4γ)
#   final 5C fit (4C + π0 mass), χ² < 40
# ---------------------------------------------------------------------------
alg_name_mode2 = "HcEtaPrimeGammapipi"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"              # at least four photons in this mode
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # final 5C kinematic fit (4C + π0 mass), loose χ² cut
  .kinematic_fit([:pi0, :gamma, :gamma, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.execute_on([psip_data, psip_incMC, exMC_mode2])

# ---------------------------------------------------------------------------
# Mode 3: η → γγ   (fully neutral final state π0 3γ)
#   final 5C fit (4C + π0 mass), χ² < 40
# ---------------------------------------------------------------------------
alg_name_mode3 = "HcEtaGG"
alg_mode3 = Algorithm.new(alg_name_mode3)
alg_mode3.set_header(["#{alg_name_mode3}Alg/#{alg_name_mode3}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:emc_relative_timing, "fully neutral final state (π0 → γγ, h_c → γ, η → γγ): " \
                                     "an additional EMC timing requirement |T - T_max| <= 500 ns " \
                                     "is applied to suppress beam-related background; this relative " \
                                     "shower-timing cut is not expressible in the DSL photon block")

sel_mode3 = Selection.new
  .select_photon {                       # fully neutral: no charged track requirement
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # final 5C kinematic fit (4C + π0 mass), loose χ² cut
  .kinematic_fit([:pi0, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)
alg_mode3.execute_on([psip_data, psip_incMC, exMC_mode3])

# ---------------------------------------------------------------------------
# Mode 4: η → π+π−π0   (charged final state π+π−π0π0γ)
#   final 6C fit (4C + 2 π0 masses), χ² < 40; competing γπ+π−π0π0 hypothesis stored as veto
# ---------------------------------------------------------------------------
alg_name_mode4 = "HcEtaPipipimpi0"
alg_mode4 = Algorithm.new(alg_name_mode4)
alg_mode4.set_header(["#{alg_name_mode4}Alg/#{alg_name_mode4}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode4 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  # tag the two π0 candidates (one from ψ(3686), one from η → π+π−π0)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # nominal 6C kinematic fit (4C + 2 π0 masses), loose χ² cut
  .kinematic_fit([:pi0, :pi0, :gamma, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
  # competing γπ+π−π0π0 hypothesis: no chi2_cut / no nominal -> stores its χ² for ROOT-level veto
  .kinematic_fit([:pi0, :pi0, :gamma, :pip, :pim]) {
    constrain_four_momentum
  }

alg_mode4.with_decay_card(decay_card_mode4).apply(sel_mode4)
alg_mode4.execute_on([psip_data, psip_incMC, exMC_mode4])

# ---------------------------------------------------------------------------
# Mode 5: h_c → γ π0 (search)   (fully neutral final state 2π0 γ)
#   final 6C fit (4C + 2 π0 masses), χ² < 20
# ---------------------------------------------------------------------------
alg_name_mode5 = "HcGammaPi0"
alg_mode5 = Algorithm.new(alg_name_mode5)
alg_mode5.set_header(["#{alg_name_mode5}Alg/#{alg_name_mode5}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:emc_relative_timing, "fully neutral final state (π0 → γγ, h_c → γ, π0 → γγ): " \
                                     "an additional EMC timing requirement |T - T_max| <= 500 ns " \
                                     "is applied to suppress beam-related background; this relative " \
                                     "shower-timing cut is not expressible in the DSL photon block")

sel_mode5 = Selection.new
  .select_photon {                       # fully neutral: no charged track requirement
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"
  }
  # tag the two π0 candidates (one from ψ(3686), one from h_c → γπ0)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # final 6C kinematic fit (4C + 2 π0 masses)
  .kinematic_fit([:pi0, :gamma, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 20
  }

alg_mode5.with_decay_card(decay_card_mode5).apply(sel_mode5)
alg_mode5.execute_on([psip_data, psip_incMC, exMC_mode5])