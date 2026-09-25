### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 phi K_S0 K_S0 PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 1M exclusive MC events for psi(3686) -> phi K_S0 K_S0
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_phiKsKs"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PhiKsKs"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy 3.686 GeV

event_selection = Selection.new
event_selection.select_track {          # charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm in the transverse plane
      nChrp     ">=3"                   # at least 3 positive tracks
      nChrn     ">=3"                   # at least 3 negative tracks
      nNet      "==0"                   # net charge zero
    }
    .pid(method: :probability) {        # kaon identification
      prob_cut 0.001                    # PID probability > 0.001
      identify :kaon, against: [:pion, :proton]  # K+ and K- (charge-conjugation shorthand)
      nkp      ">=1"                    # at least one K+
      nkm      ">=1"                    # at least one K-
    }
    .remove([:kp <= :chrgp, :km <= :chrgn])  # drop identified kaons from charged track lists
    .pid(method: :probability) {        # pion identification on the remaining tracks
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]  # pi+ and pi-
      npip     ">=2"                    # at least two pi+
      npim     ">=2"                    # at least two pi-
    }
    .secondary_vertex_fit([:pip, :pim]) {  # first K_S0 -> pi+ pi- vertex
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {  # second K_S0 -> pi+ pi- vertex
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
    # 4C kinematic fit to phi K_S0 K_S0 (phi -> K+ K-), constraining the final state to CMS 4-momentum
    .kinematic_fit([:kp, :km, :K_S0, :K_S0]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:kp, :km).within(1.01, 1.03)  # phi -> K+ K- mass window
      chi2_cut 50
    }

my_algorithm
  .note(:helix_correction,
        "helix-parameter correction applied to all charged tracks before the 4C kinematic fit; the efficiency difference with/without the correction is estimated by re-running the BOSS selection")
  .note(:ks_mass_window,
        "K_S0 candidates restricted to the mass window 0.486-0.510 GeV/c^2; sidebands 0.454-0.478 and 0.518-0.542 GeV/c^2 are used for background estimation (sidebands filled in the ROOT analysis)")
  .note(:ks_decay_length,
        "K_S0 decay length required to exceed 2 sigma and the secondary-track |Vz| < 20 cm; applied on the BOSS side after the secondary vertex fit")
  .note(:continuum_background,
        "continuum QED background estimated from off-resonance data rescaled to the psi(3686) energy by a 1/s factor")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])