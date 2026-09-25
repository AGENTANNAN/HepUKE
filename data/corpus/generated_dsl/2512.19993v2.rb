# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data, E_cm = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the signal process:
#   ψ(3686) → γ χ_c1 ,  χ_c1 → p p̄ K_S0 K_S0 ,  each K_S0 → π+π−
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 1,000,000 signal exclusive MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c1_ppbar_KS0KS0"
  config.related_dataset = psip_data          # matched to the real data taking conditions
  config.events          = 1_000_000          # 1e6 events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: "temp_for_test")

### Event selection (BOSS) ###
alg_name = "GamChiCJppbarKSKS"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # centre-of-mass energy 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                       # charged-track selection
     cos_theta 0.93                     # |cosθ| < 0.93
     Vz        10.0                     # |Vz| < 10 cm
     Vr        1.0                      # Vr < 1 cm
     nChrp     ">=3"                    # at least 3 positively charged tracks
     nChrn     ">=3"                    # at least 3 negatively charged tracks
     nNet      "==0"                    # net charge zero
  }
  .select_photon {                      # photon selection
     tdc_emc_start     0                # EMC timing window 0–14
     tdc_emc_end       14
     angle_to_track    10.0             # > 10° from the nearest charged track
     energyThreshold_b 0.025            # barrel threshold 25 MeV
     energyThreshold_e 0.050            # endcap threshold 50 MeV
     nGam              ">=1"            # at least one photon
  }
  .pid(method: :probability) {          # PID with the probability method
     prob_cut 0.001                     # PID probability > 0.001
     identify :proton, against: [:kaon, :pion]   # p / p̄ separated from K and π
     nprp ">=1"                         # at least one proton
     nprm ">=1"                         # at least one anti-proton
  }
  .remove([:prp <= :chrgp])             # protons are no longer generic charged tracks
  .remove([:prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})       # remaining tracks are taken as pions
  .secondary_vertex_fit([:pip, :pim]) {           # 1st K_S0 → π+π− secondary vertex
     build_virtual_particle(:K_S0).by_minimizing_mass_difference   # pick the π+π− pair closest to m_K0
     remove_used_particle_from_candidate_list                     # pions are consumed
  }
  .secondary_vertex_fit([:pip, :pim]) {           # 2nd K_S0 → π+π− secondary vertex
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :prp, :prm, :K_S0, :K_S0]) {   # 4C fit to γ p p̄ K_S0 K_S0
     nominal                        # nominal fit — only its corrected four-momenta are stored
     constrain_four_momentum        # 4C energy–momentum conservation against the CMS
     chi2_cut 200                   # loose BOSS value; the paper's tighter <50 is optimised in ROOT
  }

# Selection criteria that have no dedicated DSL expression, preserved for the downstream stages
my_algorithm
  .note(:ks_selection,
        "each K_S0 candidate is required to satisfy the K_S0 mass window " \
        "|M(π+π−) − m_K0| < 0.012 GeV/c² and a decay length > 2σ from the interaction point")
  .note(:background_veto,
        "Λ / Λ̄ contamination is removed by requiring |M(pπ−) − m_Λ| > 0.012 GeV/c² and " \
        "|M(p̄π+) − m_Λ̄| > 0.012 GeV/c²")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Run the selection on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])