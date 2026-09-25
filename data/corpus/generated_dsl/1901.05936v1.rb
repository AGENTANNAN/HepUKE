### Dataset preparation ###
# ψ(3770) data and inclusive MC (sample name convention: [BOSS_version]_[CMS_energy_MeV])
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the full ψ(3770) → D+ D- signal chain (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-   PHSP;
    Enddecay

    Decay D+
    1.000  K_S0  pi+  pi+  pi-   PHSP;
    Enddecay

    Decay D-
    1.000  K+  pi-  pi-   PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DplusKsPiPiPi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "DplusToKsPiPiPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770) CMS energy in GeV

event_selection = Selection.new
  .select_track {                 # charged track selection
      cos_theta 0.93              # |cosθ| < 0.93
      Vz        20.0              # |Vz| < 20 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     ">=3"             # at least 3 positive tracks
      nChrn     ">=3"             # at least 3 negative tracks
      nNet      "==0"             # net charge zero
  }
  .pid(method: :probability) {    # identify the K+ of the D- tag first
      prob_cut 0.001
      identify :kaon, against: [:pion]
      nkp      ">=1"              # at least one kaon
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])   # remove identified kaons from the charged lists
  .pid(method: :probability) {    # then identify the remaining tracks as pions
      prob_cut 0.001
      identify :pion, against: [:kaon]
  }
  .secondary_vertex_fit([:pip, :pim]) {     # build K_S0 from π+π- pairs
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # 6C kinematic fit: 4C + m(D+) + m(K_S0)
  .kinematic_fit([:K_S0, :pip, :pip, :pim, :kp, :pim, :pim]) {
      nominal
      invariant_mass_of(:K_S0).within(0.4676, 0.5276)                            # K_S0 mass window
      constrain_four_momentum                                                    # 4C energy-momentum conservation
      invariant_mass_of(:K_S0, :pip, :pip, :pim).constrain_to_nominal_mass_of(:Dplus)   # D+ mass constraint
      invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)               # K_S0 mass constraint
      chi2_cut 100
  }

# BOSS-side procedures that have no dedicated DSL construct
my_algorithm
  .note(:delta_e_window, "ΔE window applied on both D candidates before the 6C fit: " \
                         "tag D- ∈ [-0.027, 0.025] GeV, signal D+ ∈ [-0.033, 0.030] GeV")
  .note(:mbc_window, "beam-constrained mass M_BC of both D candidates required in " \
                     "[1.8628, 1.8788] GeV/c² before the 6C fit")
  .note(:background_veto, "events containing an additional K_S0 candidate are vetoed: an extra " \
                          "π+π- pair whose invariant mass is within 30 MeV/c² of the K_S0 nominal " \
                          "mass and whose decay-length significance exceeds 2σ")
  .note(:prefit_vertex_fit, "a preliminary vertex fit on all charged tracks with χ² < 100 is " \
                            "performed before the 6C kinematic fit to improve the track parameters")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute the algorithm on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([data_3773, incMC_3773, exMC_signal])