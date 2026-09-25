# Core DSL classes and dependencies are loaded automatically at execution time.

### Dataset preparation ###
# 31 scan energy points spanning 3.650 - 4.914 GeV (real data + matching inclusive MC).
scan_samples = [
  "709_3650", "709_3686",
  "703_4009", "705_4130", "705_4160", "703_4180", "703_4190",
  "703_4200", "703_4210", "703_4220", "703_4230", "703_4237",
  "703_4246", "703_4260", "703_4270", "703_4280", "703_4360",
  "703_4420", "703_4600", "706_4610", "706_4620", "706_4640",
  "706_4660", "706_4680", "706_4700", "707_4740", "707_4750",
  "707_4780", "707_4840", "707_4914", "707_4946"
]
scan_data  = scan_samples.map { |s| DatasetManager.real_data.find(s) }       # real data at every scan point
scan_incMC = scan_samples.map { |s| DatasetManager.inclusive_mc.find(s) }   # corresponding inclusive MC

# Decay card for e+e- -> K- anti-Xi- Lambda; ISR / vacuum polarisation modelled by the ConExc
# generator (no KKMC + psi(4260) default, no explicit "Particle vpho" - injected per energy point).
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 K- anti-Xi- Lambda ConExc;
    Enddecay

    Decay anti-Xi-
    1.000 anti-Lambda pi+ PHSP;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC generated for every energy point (one MC per dataset).
exMCs = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_KXiLambda"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: "temp_KXiLambda") }

### Event selection (BOSS) ###
alg_name = "KXiLambda"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.650]})            # nominal beam energy of the scan
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                                                # charged track selection
        cos_theta 0.93          # |cos(theta)| < 0.93
        nChrp     ">=2"         # at least two positive tracks
        nChrn     ">=2"         # at least two negative tracks
    }
    .pid(method: :probability) {                                   # PID by the probability method
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]                  # proton vs kaon/pion separation
        identify :kaon,   against: [:pion]                         # kaon vs pion separation
        identify :pion,   against: [:kaon, :proton]                # pion list
        nprm ">=1"          # at least one anti-proton
        nkm  ">=1"          # at least one K-
        npip ">=2"          # at least two pi+
    }
    .secondary_vertex_fit([:prm, :pip]) {                          # reconstruct anti-Lambda -> anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list                   # do not reuse the fitted tracks
    }
    .partial_miss([3]) {                                           # partial reconstruction: Lambda (recID 3) is missing
        best_combination_by_mass :Lambda_bar,    1.11568           # best (anti-Lambda, pi+) candidate by combined mass difference
        best_combination_by_mass :Xi_bar_minus,  1.32171           # form anti-Xi- -> anti-Lambda pi+
        require_recoil_mass 1.0, 1.3                               # recoil mass against K- anti-Xi- in 1.0-1.3 GeV/c^2
    }

my_Algorithm
    .note(:mass_window, "anti-Lambda -> anti-p pi+ required within +-8 MeV/c^2 of the nominal Lambda mass; "
                      + "anti-Xi- -> anti-Lambda pi+ within +-6 MeV/c^2 of the nominal anti-Xi- mass; the "
                      + "best (anti-Lambda, anti-Xi-) combination is chosen by minimising the combined mass difference")
    .note(:decay_length, "the reconstructed anti-Lambda secondary vertex is required to have a positive decay length")
    .note(:isr_correction, "ISR correction is obtained by iterative QED reweighting of the signal MC (ConExc generator)")
    .note(:sigma0_mode, "the K- anti-Xi- Sigma0 (Sigma0 -> Lambda gamma) mode is implicitly selected by the same "
                      + "1.0-1.3 GeV/c^2 recoil-mass window, which covers both the Lambda and Sigma0 peaks")
    .note(:charge_conjugation, "the charge-conjugate channel e+e- -> K+ Xi- anti-Lambda is implicitly included")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on all 31 scan points together with the generated exclusive signal MC.
root_files = my_Algorithm.execute_on(scan_data + scan_incMC + exMCs)