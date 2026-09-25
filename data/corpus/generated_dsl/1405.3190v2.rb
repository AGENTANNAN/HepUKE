# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data (~2.92 fb⁻¹)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # inclusive MC at 3.773 GeV

# --- Decay cards for the three radiative signal modes ---
# ψ(3770) → γ η_c, η_c → K_S0 K+ π− (charge conjugate implied), K_S0 → π+ π−
decay_card_etac = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ψ(3770) → γ η_c(2S), η_c(2S) → K_S0 K+ π−
decay_card_etac2s = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ψ(3770) → γ χ_c1 (cross-check), χ_c1 → K_S0 K+ π−
decay_card_chic1 = <<~DECAYCARD
    Decay psi(3770)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC, 100k events for each of the three signal modes ---
exMC_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_gamma_etac_ks0_kpi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_etac
  config.cross_section   = :default
end

exMC_etac2s = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_gamma_etac2s_ks0_kpi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_etac2s
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_gamma_chic1_ks0_kpi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# The three modes share identical final state (γ K_S0 K+π−) and selection -> one Algorithm.
alg_name = "RadiativeKpi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770) CMS energy in GeV

# Common selection chain shared by all three signal modes.
event_selection = Selection.new
    .select_track {                  # charged-track selection
        cos_theta 0.93               # |cosθ| < 0.93
        Vz        10.0               # |Vz| < 10 cm
        Vr        10.0               # Vr < 1 cm (= 10 mm)
        nChrp     ">=2"              # ≥ 2 positive tracks
        nChrn     ">=2"              # ≥ 2 negative tracks  (≥ 4 charged tracks in total)
        nNet      "==0"              # net charge zero
    }
    .select_photon {                 # photon selection
        tdc_emc_start     0
        tdc_emc_end       14         # EMC timing within 700 ns of the collision
        angle_to_track    20.0       # > 20° from any charged track
        energyThreshold_b 0.025      # E > 25 MeV in the barrel
        energyThreshold_e 0.050      # E > 50 MeV in the endcap
        nGam              ">=1"      # at least one photon
    }
    .assign({:chrgp => :pip, :chrgn => :pim})    # treat the K_S0 daughters as pions
    .secondary_vertex_fit([:pip, :pim]) {        # reconstruct K_S0 → π+π−
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .pid(method: :chi2_sum) {                    # χ²-sum PID for the remaining K and π
        chi_min_cut 4                            # χ²_PID < 4
        identify :kaon, :pion
    }
    .kinematic_fit([:gamma, :K_S0, :kp, :pim]) { # nominal 4C fit to γ K_S0 K+π− (and c.c.)
        nominal
        constrain_four_momentum
        chi2_cut 200                             # loose in BOSS; paper applies χ²_4C < 20 in ROOT
    }
    .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pim]) {  # 2nd 4C fit (γγ hypothesis): stores χ² for π0 veto
        constrain_four_momentum
    }

# Inexpressible BOSS-side procedures captured as notes.
my_algorithm
    .note(:ks0_selection,
          "K_S0 candidates required to satisfy |M(π+π−) − M(K_S0)| < 10 MeV/c² " \
          "(window 0.4876–0.5076 GeV/c²) and decay length > 2× vertex resolution; " \
          "these windows are not expressible in the vertex-fit block")
    .note(:background_veto,
          "D0 veto applied: events with |M(K±π∓) − M(D0)| < 3σ rejected to suppress " \
          "ψ(3770) → D0 D̄0 background")
    .note(:helix_correction,
          "charged-track helix parameters corrected for data/MC differences using the " \
          "J/ψ → φ f0(980) control sample")
    .with_decay_card(decay_card_etac)
    .apply(event_selection)

# Execute on real data, inclusive MC and the three signal exclusive MC samples
root_files = my_algorithm.execute_on([data_3773, incMC_3773, exMC_etac, exMC_etac2s, exMC_chic1])