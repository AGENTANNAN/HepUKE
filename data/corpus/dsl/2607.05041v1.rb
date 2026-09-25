# Dataset preparation
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> gamma chi_cJ, chi_cJ -> p K- anti-Lambda eta (+c.c.)
# One card per chi_cJ. All three share the same final state and selection.
decay_card_chic0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0                               PHSP;
  Enddecay

  Decay chi_c0
  1.0000 p+ K- anti-Lambda0 eta                     PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                                PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                                PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1                               PHSP;
  Enddecay

  Decay chi_c1
  1.0000 p+ K- anti-Lambda0 eta                     PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                                PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                                PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2                               PHSP;
  Enddecay

  Decay chi_c2
  1.0000 p+ K- anti-Lambda0 eta                     PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                                PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                                PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC samples for each chi_cJ signal mode
exMC_chic0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "chic0_pKLambdaEta"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card_chic0
  c.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "chic1_pKLambdaEta"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card_chic1
  c.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "chic2_pKLambdaEta"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card_chic2
  c.cross_section   = :default
end

# Event selection: chi_cJ -> p K- anti-Lambda eta with anti-Lambda -> anti-p pi+, eta -> gamma gamma
alg_name    = "ChicJpKLambdaEta"
algorithm   = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93 for MDC tracks
    Vz        10.0          # |Vz| < 10 cm for tracks from chi_cJ decay
    Vr        1.0           # |Vxy| < 1 cm for tracks from chi_cJ decay
    nChrp     ">=2"         # at least: p, pi+ (from anti-Lambda)
    nChrn     ">=2"         # at least: K-, anti-p (from anti-Lambda)
    nNet      "==0"         # net charge = 0
  }
  .select_photon {
    tdc_emc_start     0     # EMC time window start (>= 0 ns)
    tdc_emc_end       14    # EMC time window end (<= 700 ns = 14 * 50 ns)
    angle_to_track    10.0  # reject photons within 10 deg of any charged track
    energyThreshold_b 0.025 # E > 25 MeV in barrel (|cos theta| < 0.80)
    energyThreshold_e 0.050 # E > 50 MeV in endcap (0.86 < |cos theta| < 0.92)
    nGam              ">=3" # at least 3 good photon candidates
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]  # L(p) > L(K) and L(p) > L(pi)
    identify :kaon,   against: [:proton, :pion] # L(K) > L(p) and L(K) > L(pi)
    nprp ">=1"
    nkm  ">=1"
  }
  .remove([:prp <= :chrgp])
  .remove([:km  <= :chrgn])
  .select_isolated_photon {
    angle_to_prm_track 20.0 # reject photons within 20 deg of anti-p track
    nGam ">=3"
  }
  # Reconstruct anti-Lambda -> anti-p pi+ via secondary vertex fit.
  # Use the remaining charged tracks (protons and kaons already removed).
  .assign({:chrgp => :pip, :chrgn => :prm})
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit: e+e- -> p K- anti-Lambda gamma gamma gamma
  .kinematic_fit([:prp, :km, :Lambda_bar, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200            # loose BOSS-side cut; tight chi2_4C < 35 applied in ROOT
  }

# Non-DSL selections captured as notes for downstream systematics:
algorithm
  .note(:lambda_secondary_vertex_chi2, "anti-Lambda secondary-vertex fit chi2_vtx < 200 required")
  .note(:lambda_decay_length,          "decay length significance L/sigma_L > 2 required for anti-Lambda")
  .note(:lambda_mass_window,           "|M(anti-p pi+) - m(Lambda)| < 0.015 GeV/c^2 for anti-Lambda signal region")
  .note(:eta_pair_selection,           "photons ordered by energy; eta from gamma1-gamma2 or gamma1-gamma3 pairs only; reject events with both pairs in eta window")
  .note(:eta_mass_window,              "eta signal window 0.525 < M(gg) < 0.570 GeV/c^2")
  .note(:pi0_veto,                     "veto |M(g1g2)-m(pi0)|<0.013, |M(g1g3)-m(pi0)|<0.018, |M(g2g3)-m(pi0)|<0.009 GeV/c^2")
  .note(:jpsi_veto,                    "veto |RM(eta) - m(J/psi)| < 0.007 GeV/c^2 (recoil mass of eta) to suppress psi(3686) -> eta J/psi")
  .note(:sigma0_veto,                  "veto |M(anti-Lambda gamma_unused) - m(Sigma0)| < 0.018 GeV/c^2 to suppress anti-Sigma0 -> anti-Lambda gamma")
  .note(:chi2_4c_tight_in_root,        "final chi2_4C < 35 requirement applied in ROOT-side analysis after FOM optimization")

algorithm.with_decay_card(decay_card_chic0).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])
