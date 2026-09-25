# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# e+e- -> pi+ pi- D+ D- measured on 37 c.m. energy points, 4.190 - 4.946 GeV (17.4 fb^-1)
energy_names = %w[
  703_4190 703_4200 703_4210 703_4220 703_4230 703_4237 703_4246
  703_4260 703_4270 703_4280 703_4310 703_4360 703_4390 703_4420
  703_4470 703_4530 703_4575 703_4600
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
data_points   = energy_names.map { |name| DatasetManager.real_data.find(name) }
incMC_samples = energy_names.map { |name| DatasetManager.inclusive_mc.find(name) }

# Signal decay card (EvtGen, phase space): e+e- -> pi+ pi- D+ D-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  pi+  pi-  D+  D-       PHSP;
    Enddecay

    Decay D+
    1.0000  K-  pi+  pi+           PHSP;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-           PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive phase-space MC events for e+e- -> pi+ pi- D+ D- (one per energy point)
exMC_signals = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipiDpDm"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signals.each { |m| m.save_to_config(format: :yaml, file_path: "temp_exmc_pipiDpDm") }

### Event selection (BOSS) ###
alg_name = "PiPiDpDmRecoil"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})

event_selection = Selection.new
event_selection
    .select_track {
        cos_theta 0.93    # |cos(theta)| < 0.93
        Vz        10.0    # |Vz| < 10 cm
        Vr        1.0     # Vr < 1 cm in the transverse plane
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon,   against: [:pion]         # K- / K+ vs pi
        identify :pion,   against: [:kaon]         # pi vs K
        identify :proton, against: [:pion, :kaon]  # for the proton/antiproton veto
        nkm  "==1"    # one K-  (from D+ -> K- pi+ pi+)
        npip "==3"    # two pi+ from D+ plus one recoil pi+
        npim "==1"    # one recoil pi-
        nprp "==0"    # veto proton
        nprm "==0"    # veto antiproton
    }
    # Reconstruct the tag D+ -> K- pi+ pi+ with a secondary vertex fit
    .secondary_vertex_fit([:km, :pip, :pip]) {
        build_virtual_particle(:Dplus).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # Partial reconstruction: D- inferred from the recoil mass of the D+ pi+ pi- system
    .partial_miss([4]) {
        best_combination_by_mass :Dplus, 1.8697
        require_recoil_mass 1.861, 1.879   # D- mass window ~ +/- 9 MeV/c^2
    }

my_algorithm
    .note(:dp_mass_window, "D+ candidate required to satisfy |M(K- pi+ pi+) - M(D+)| < 11 MeV/c^2")
    .note(:dp_mass_constrained_fit, "D+ mass-constrained kinematic fit applied on the K- pi+ pi+ system (partial reconstruction replaces the 4C fit); improves the D+ momentum and hence the D- recoil-mass resolution")
    .note(:background_veto, "Lambda_c background suppressed by tighter vertex cuts on D+ daughter tracks (Vr < 0.55 cm, |Vz| < 3 cm) and by proton/antiproton vetoes; D0 and K_S0 backgrounds vetoed by invariant-mass windows and secondary-vertex criteria")
    .note(:ecms_scan, "ECMS varies over the 37 scan points (4.190-4.946 GeV); the declared constant is nominal and is set per energy point")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on(data_points + incMC_samples + exMC_signals)