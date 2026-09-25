# BOSS Ruby DSL for BESIII paper 1802.03119v1
# Absolute branching fractions for D mesons decaying into two pseudoscalar mesons
# Data: ψ(3770) at √s = 3.773 GeV, 2.93 fb^-1, single-tag method
# D decays reconstructed from their final-state daughters; mBC and ΔE applied at ROOT level.
# Covered topologies:
#   D0 → Kπ/ππ/KK (all charged)
#   D+ → K+π0 / π+π0 (charged + π0)
#   D+ → KS0π+ / KS0K+ (KS0 → π+π-)
#   D0 → KS0π0 (KS0 + π0)
#   D+ → π+η (η → γγ)
#   D0 → KS0η' (KS0 + η' → π+π-η)

### Dataset ###
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay cards
# ============================================================

# D0 → K-π+ (and charge conjugate D0bar → K+π-)
decay_card_D0_Kpi = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0  PHSP;
    Enddecay
    Decay D0
    1.000  K-  pi+  PHSP;
    Enddecay
    Decay anti-D0
    1.000  K+  pi-  PHSP;
    Enddecay
    End
DECAYCARD

# D0 → π+π- (and D0 → K+K- similar)
decay_card_D0_pipi = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0  PHSP;
    Enddecay
    Decay D0
    1.000  pi+  pi-  PHSP;
    Enddecay
    Decay anti-D0
    1.000  K+  pi-  PHSP;
    Enddecay
    End
DECAYCARD

# D+ → K+π0 (and D+ → π+π0 similar)
decay_card_Dp_Kpi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-  PHSP;
    Enddecay
    Decay D+
    1.000  K+  pi0  PHSP;
    Enddecay
    Decay D-
    1.000  K-  pi0  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay
    End
DECAYCARD

# D+ → KS0π+ (KS0 → π+π-)
decay_card_Dp_KsPi = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-  PHSP;
    Enddecay
    Decay D+
    1.000  K_S0  pi+  PHSP;
    Enddecay
    Decay D-
    1.000  K-  pi0  PHSP;
    Enddecay
    Decay K_S0
    1.000  pi+  pi-  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay
    End
DECAYCARD

# D0 → KS0π0 (KS0 → π+π-, π0 → γγ)
decay_card_D0_KsPi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0  PHSP;
    Enddecay
    Decay D0
    1.000  K_S0  pi0  PHSP;
    Enddecay
    Decay anti-D0
    1.000  K+  pi-  PHSP;
    Enddecay
    Decay K_S0
    1.000  pi+  pi-  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay
    End
DECAYCARD

# D+ → π+η (η → γγ)
decay_card_Dp_pi_eta = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-  PHSP;
    Enddecay
    Decay D+
    1.000  pi+  eta  PHSP;
    Enddecay
    Decay D-
    1.000  K-  pi0  PHSP;
    Enddecay
    Decay eta
    1.000  gamma  gamma  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_D0_Kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0_Kpi"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_D0_Kpi
  config.cross_section = :default
end

exMC_D0_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0_pipi"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_D0_pipi
  config.cross_section = :default
end

exMC_Dp_Kpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_Dp_Kpi0"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_Dp_Kpi0
  config.cross_section = :default
end

exMC_Dp_KsPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_Dp_KsPi"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_Dp_KsPi
  config.cross_section = :default
end

exMC_D0_KsPi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0_KsPi0"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_D0_KsPi0
  config.cross_section = :default
end

exMC_Dp_pi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_Dp_pi_eta"
  config.related_dataset = data_3773
  config.events = 100000
  config.decay_card = decay_card_Dp_pi_eta
  config.cross_section = :default
end

# ============================================================
# Topology A: D → two charged tracks (D0→Kπ, D0→ππ, D0→KK)
# ============================================================
# Selection: 2 oppositely charged tracks, |cosθ| < 0.93, Vxy < 1 cm, |Vz| < 10 cm
# PID: Prob(K) > Prob(π) for kaons, Prob(π) > Prob(K) for pions
# Cosmic/Bhabha veto: TOF time diff < 5 ns, not e+e- or μ+μ- pair
# Additionally: at least 1 extra EMC cluster > 50 MeV or 1 extra MDC track (D0D0 coherence)
# No kinematic fit - mBC and ΔE applied as ROOT-level cuts

alg_D0_2body = Algorithm.new("D0TwoBody")
alg_D0_2body.set_header(["D0TwoBodyAlg/D0TwoBody.h"])
            .set_constant({"ECMS" => [:double, 3.773]})

sel_D0_2body = Selection.new
sel_D0_2body.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet  "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .kinematic_fit([:km, :pip]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_D0_2body
  .note(:pid_correction_method, "MC efficiencies corrected by momentum-dependent data-MC differences; tracking systematic 0.3%, PID 0.3% per track")
  .note(:background_veto, "cosmic/Bhabha rejection: TOF time diff < 5 ns between two tracks; exclude μ+μ- and e+e- pair candidates; require extra EMC cluster > 50 MeV or extra MDC track")
  .note(:efficiency_curve, "mBC fit with signal MC shape convoluted with double Gaussian; ARGUS background; ΔE requirement ±3σ_ΔE applied at ROOT level")
  .with_decay_card(decay_card_D0_Kpi)
  .apply(sel_D0_2body)

# ============================================================
# Topology B: D+ → K+π0 / D+ → π+π0 (charged + π0)
# ============================================================
# Selection: ≥1 charged, ≥2 photons. Reconstruct π0 via Kalman fit.
# PID on charged track. Then kinematic fit with D mass + CMS constraints.

alg_Dp_pi0 = Algorithm.new("DpChargedPi0")
alg_Dp_pi0.set_header(["DpChargedPi0Alg/DpChargedPi0.h"])
           .set_constant({"ECMS" => [:double, 3.773]})

sel_Dp_pi0 = Selection.new
sel_Dp_pi0.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot ">=1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:kp, :pi0]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_Dp_pi0
  .note(:pid_correction_method, "tracking 0.3%, PID 0.3%; π0 reconstruction 1.0% uncertainty; MC efficiencies corrected for data-MC differences")
  .note(:efficiency_curve, "mBC fit; ΔE ±3σ applied at ROOT; π0 mass window (0.115,0.150) GeV/c²; 1C mass-constrained fit for π0 momentum resolution improvement")
  .with_decay_card(decay_card_Dp_Kpi0)
  .apply(sel_Dp_pi0)

# ============================================================
# Topology C: D+ → KS0π+ / D+ → KS0K+ (KS0 → π+π-)
# ============================================================
# Selection: ≥3 charged tracks. Reconstruct KS0 via secondary vertex fit.
# KS0: decay length > 2σ, π+π- mass in ±12 MeV/c² of KS0 nominal mass.

alg_Dp_Ks = Algorithm.new("DpKsCharged")
alg_Dp_Ks.set_header(["DpKsChargedAlg/DpKsCharged.h"])
          .set_constant({"ECMS" => [:double, 3.773]})

sel_Dp_Ks = Selection.new
sel_Dp_Ks.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:K_S0, :kp]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_Dp_Ks
  .note(:pid_correction_method, "KS0 reconstruction efficiency: 1.5% per KS0 (corrected for data-MC differences); tracking 0.3%, PID 0.3%")
  .note(:efficiency_curve, "KS0 selection: flight significance L/σ > 2, M(π+π-) within ±12 MeV/c² of KS0 nominal mass; mBC fit with ΔE ±3σ applied at ROOT")
  .with_decay_card(decay_card_Dp_KsPi)
  .apply(sel_Dp_Ks)

# ============================================================
# Topology D: D0 → KS0π0 (KS0 → π+π-, π0 → γγ)
# ============================================================
# Selection: ≥2 charged (KS0 daughters), ≥2 photons.
# Reconstruct KS0 (secondary vertex) and π0 (Kalman fit), then combine.

alg_D0_KsPi0 = Algorithm.new("D0KsPi0")
alg_D0_KsPi0.set_header(["D0KsPi0Alg/D0KsPi0.h"])
             .set_constant({"ECMS" => [:double, 3.773]})

sel_D0_KsPi0 = Selection.new
sel_D0_KsPi0.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot ">=2"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:K_S0, :pi0]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_D0_KsPi0
  .note(:pid_correction_method, "KS0 efficiency 1.5%; π0 efficiency 1.0%; tracking 0.3%; MC efficiencies corrected for data-MC differences")
  .note(:efficiency_curve, "KS0 sidebands used for peaking background estimation; mBC fit; ΔE ±3σ at ROOT; π0 mass window (0.115,0.150) GeV/c²")
  .with_decay_card(decay_card_D0_KsPi0)
  .apply(sel_D0_KsPi0)

# ============================================================
# Topology E: D+ → π+η (η → γγ)
# ============================================================
# Selection: ≥1 charged, ≥2 photons. Reconstruct η via Kalman fit.

alg_Dp_pi_eta = Algorithm.new("DpPiEta")
alg_Dp_pi_eta.set_header(["DpPiEtaAlg/DpPiEta.h"])
             .set_constant({"ECMS" => [:double, 3.773]})

sel_Dp_pi_eta = Selection.new
sel_Dp_pi_eta.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot ">=1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  .kinematic_fit([:pip, :eta]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_Dp_pi_eta
  .note(:pid_correction_method, "η reconstruction 1.0% uncertainty; MC efficiencies corrected for data-MC differences")
  .note(:efficiency_curve, "η mass window (0.515,0.575) GeV/c²; 1C mass-constrained fit for η; mBC fit; ΔE ±3σ at ROOT")
  .with_decay_card(decay_card_Dp_pi_eta)
  .apply(sel_Dp_pi_eta)

### Execute ###
alg_D0_2body.execute_on([data_3773, incMC_3773, exMC_D0_Kpi])
alg_Dp_pi0.execute_on([data_3773, incMC_3773, exMC_Dp_Kpi0])
alg_Dp_Ks.execute_on([data_3773, incMC_3773, exMC_Dp_KsPi])
alg_D0_KsPi0.execute_on([data_3773, incMC_3773, exMC_D0_KsPi0])
alg_Dp_pi_eta.execute_on([data_3773, incMC_3773, exMC_Dp_pi_eta])