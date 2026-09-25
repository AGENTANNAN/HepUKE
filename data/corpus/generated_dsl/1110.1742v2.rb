# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")        # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # Corresponding inclusive MC

# Decay card: ψ(2S) → γ χ_c2 , χ_c2 → π+π−  (phase space)
decay_card_pipi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay
    Decay chi_c2
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: ψ(2S) → γ χ_c2 , χ_c2 → K+K−  (phase space)
decay_card_KK = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay
    Decay chi_c2
    1.000 K+ K- PHSP;
    Enddecay
    End
DECAYCARD

# 200k-event phase-space exclusive MC for each mode (used for normalisation / PWA)
exMC_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gchic2_pipi"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_pipi
  config.cross_section  = :default
end

exMC_KK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gchic2_KK"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_KK
  config.cross_section  = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------
# Mode I : ψ(2S) → γ χ_c2 , χ_c2 → π+π−
# ------------------------------------------------------------------
alg_name_pipi = "PsipToGChiC2PiPi"
alg_pipi = Algorithm.new(alg_name_pipi)
alg_pipi.set_header(["#{alg_name_pipi}Alg/#{alg_name_pipi}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        # Inexpressible BOSS-side background suppression
        .note(:background_veto, "electron background suppressed: each track EMC energy < 1.4 GeV, " \
                                "and dE/dx within 3 sigma of the expected value (tightened to 2 sigma " \
                                "in the EMC-insensitive region 0.81 < |cos(theta)| < 0.86)")
        .note(:muon_veto, "at least one track required to have EMC energy > 0.34 GeV to reject " \
                          "muon-pair backgrounds")
        .note(:radiative_photon, "the highest-energy selected photon is taken as the radiative gamma")

sel_pipi = Selection.new
  .select_track {                       # exactly one π+ and one π−, net charge zero
    cos_theta 0.93                      # |cosθ| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                      # radiative photon: E>25 MeV barrel / E>50 MeV endcap, gap excluded automatically
    tdc_emc_start   0
    tdc_emc_end     14
    angle_to_track  20.0                # > 20° from any charged track
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {          # ππ mode: one π+ and one π− against K and p
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  # Nominal 4C fit under the γ π+π− hypothesis (loose χ² cut)
  .kinematic_fit([:gamma, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Competing γ K+K− hypothesis: store its χ² for later arbitration (no cut, not nominal)
  .assign({:pip => :kp, :pim => :km})
  .kinematic_fit([:gamma, :kp, :km]) {
    constrain_four_momentum
  }

alg_pipi.with_decay_card(decay_card_pipi).apply(sel_pipi)

# ------------------------------------------------------------------
# Mode II : ψ(2S) → γ χ_c2 , χ_c2 → K+K−
# ------------------------------------------------------------------
alg_name_KK = "PsipToGChiC2KK"
alg_KK = Algorithm.new(alg_name_KK)
alg_KK.set_header(["#{alg_name_KK}Alg/#{alg_name_KK}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:background_veto, "electron background suppressed: each track EMC energy < 1.4 GeV, " \
                              "and dE/dx within 3 sigma of the expected value")
      .note(:radiative_photon, "the highest-energy selected photon is taken as the radiative gamma")

sel_KK = Selection.new
  .select_track {                       # exactly one K+ and one K−, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start   0
    tdc_emc_end     14
    angle_to_track  20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {          # KK mode: one K+ and one K− , Prob(K)>Prob(π) and Prob(K)>0.001
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  }
  # Nominal 4C fit under the γ K+K− hypothesis (loose χ² cut)
  .kinematic_fit([:gamma, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Competing γ π+π− hypothesis: store its χ² for later arbitration (no cut, not nominal)
  .assign({:kp => :pip, :km => :pim})
  .kinematic_fit([:gamma, :pip, :pim]) {
    constrain_four_momentum
  }

alg_KK.with_decay_card(decay_card_KK).apply(sel_KK)

### Execute on datasets ###
# NOTE: the final tight selection (χ² < 60 and the smaller χ² of the two hypotheses)
# is applied in the ROOT analysis, not here.
alg_pipi.execute_on([psip_data, psip_incMC, exMC_pipi])
alg_KK.execute_on([psip_data, psip_incMC, exMC_KK])