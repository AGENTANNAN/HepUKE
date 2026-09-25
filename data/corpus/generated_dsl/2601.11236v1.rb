# ============================================================================
# ψ(3686) → π0 h_c,  h_c → γ η_c,  η_c → γγ,  π0 → γγ
# BOSS part only: dataset preparation + event selection up to the 4C fit
# ============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# Decay card for the signal chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c      PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c  PHSP;
    Enddecay

    Decay eta_c
    1.0000 gamma gamma  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events for this chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_gamma_etac"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name     = "Pi0HcEtac"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new

# Charged-track selection: tracks are needed for the π+π−, π0π0 and π+π−π0 background vetoes
event_selection.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  }
  # ---------------- Photon selection ----------------
  .select_photon {
    tdc_emc_start     0        # EMC timing window start (within 700 ns)
    tdc_emc_end       14       # EMC timing window end (within 700 ns)
    angle_to_track    10.0     # at least 10° from any charged track
    energyThreshold_b 0.025    # E > 25 MeV in the barrel
    energyThreshold_e 0.050    # E > 50 MeV in the endcap
    nGam              ">=4"    # at least four photons in the event
  }
  # ---------------- Pion identification (for the charged-track vetoes) ----------------
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # π+ and π− (charge-conjugation shorthand)
  }
  # ---------------- Tag-side π0: 1C mass-constrained Kalman fit from two photons ----------------
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # constrain M(γγ) to m(π0)
    chi2_cut 5     # require the 1C fit χ² < 5
    npi0     ">=1" # at least one π0 candidate
  }
  # ---------------- Nominal 4C kinematic fit of π0 + 3γ (E1 photon + η_c → γγ) ----------------
  .kinematic_fit([:pi0, :gamma, :gamma, :gamma]) {
    nominal                                                              # nominal fit: corrected 4-momenta used
    constrain_four_momentum                                              # constrain to the initial e+e− four-momentum
    invariant_mass_of(:gamma, :gamma).within(2.8891, 3.0668)             # M(η_c → γγ) window (GeV/c²)
    chi2_cut 70                                                          # require χ²(4C) < 70
  }

# --- BOSS-side procedures that cannot be expressed in the current DSL ---
my_algorithm
  .note(:efficiency_curve,
        "tag-side π0 is built only from two BARREL photons with E > 40 MeV, pre-selected in " \
        "0.120 < M(γγ) < 0.145 GeV/c² before the 1C mass constraint; when several π0 " \
        "candidates survive, the smallest-χ²(1C) combination is kept; the π0 recoil mass is " \
        "required in the h_c window [3.48, 3.57] GeV/c²")
  .note(:efficiency_curve,
        "E1 photon (from h_c → γ η_c) is chosen among the photons not used by the tag-side " \
        "π0: any photon that pairs with another photon to form a π0 candidate " \
        "(M(γγ) in [0.115, 0.150] GeV/c² with 1C χ² < 200) is rejected; the E1 photon energy " \
        "measured in the π0-recoiling system is required in [0.46, 0.58] GeV; one candidate " \
        "is picked at random when several remain")
  .note(:background_veto,
        "π+π− recoil-mass veto: events with M_recoil(π+π−) within ±7 MeV of m(J/ψ) are rejected")
  .note(:background_veto,
        "π0π0 recoil-mass veto: events with M_recoil(π0π0) in [m(J/ψ)−15, m(J/ψ)+25] MeV are rejected")
  .note(:background_veto,
        "π+π−π0 invariant-mass veto: events with M(π+π−π0) within ±12 MeV of m(η) are rejected")
  .note(:background_veto,
        "signal-side veto: reject events where the E1γ is paired with an η_c γ in the " \
        "[0.50, 0.56] GeV or [0.88, 0.98] GeV window")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])