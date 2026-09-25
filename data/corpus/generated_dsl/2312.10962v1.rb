# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data: representative points of the 3.51-4.95 GeV energy scan
data_3773 = DatasetManager.real_data.find("712_3773")   # 3.773 GeV (psi(3770))
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.180 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.230 GeV
data_4914 = DatasetManager.real_data.find("707_4914")   # 4.914 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # 4.946 GeV
scan_data = [data_3773, data_4180, data_4230, data_4914, data_4946]

# Matching inclusive MC samples
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")
scan_incMC = [incMC_3773, incMC_4180, incMC_4230, incMC_4914, incMC_4946]

# ConExc decay card: continuum e+e- -> K_S0 K_L0 (VSS: vector -> scalar + scalar),
# ISR modelled to second order. The DSL auto-detects the literal ConExc token and
# switches to the no-KKMC template, injecting "Particle vpho <ECMS> 0.0" per energy
# point - therefore NO Particle vpho line is written for this multi-energy scan.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 8;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC generated at every scan point (shared card / cross section / events)
exMCs = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "kskl_scan_exclusive_mc"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KSKL"
kskl_alg = Algorithm.new(alg_name)
kskl_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})   # nominal reference; per-point beam energy read at run time
        .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {
    nChrp "==1"     # exactly one positive track
    nChrn "==1"     # exactly one negative track
    nNet  "==0"     # net charge zero (one pi+ and one pi-)
    # no |cos(theta)|, |Vz| or |Vr| cuts encoded
  }
  .select_photon {
    tdc_emc_start      0        # EMC timing window 0-700 ns (start)
    tdc_emc_end        14       # 14 x 50 ns = 700 ns (end)
    energyThreshold_b  0.025    # E_gamma > 25 MeV (barrel)
    energyThreshold_e  0.025    # E_gamma > 25 MeV (endcap)
    angle_to_track     20.0     # > 20 degrees from the nearest charged track
  }
  .pid(method: :probability) {
    prob_cut 0.001                                              # PID probability > 0.001
    identify :pion, against: [:kaon, :proton, :electron, :muon]  # pi+/- over K, p, e, mu
    npip "==1"
    npim "==1"
  }
  .secondary_vertex_fit([:pip, :pim]) {           # K_S0 -> pi+ pi- secondary vertex
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
  }
  .kinematic_fit([:pip, :pim]) {                  # 1C fit: m(pi+pi-) = m(K_S0)
    nominal
    invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
    chi2_cut 12
  }

kskl_alg
  .note(:ks_selection, "K_S0 secondary-vertex-fit selection: m(pi+pi-) in [0.478, 0.518] GeV/c^2, K_S0 decay length > 2 cm and vertex-fit chi^2 < 15 (these cut values are not expressible inside the secondary_vertex_fit block)")
  .note(:pid_correction_method, "electron rejection beyond the PID likelihood requirement L(pi) > L(e): each track must satisfy E/p < 1.2 GeV/(GeV/c) (E/p is not an expressible DSL cut)")
  .note(:background_veto, "reject the event if any photon pair has M(gamma gamma) in [0.123, 0.144] GeV/c^2 (pi0 -> gamma gamma veto); require the total neutral (EMC) energy outside the 20-degree cone opposite the K_S0 direction to be < 0.2 GeV to suppress e+e- -> K*0(892) K0 + c.c.")
  .note(:kl_reconstruction, "K_L0 is not explicitly reconstructed; its momentum is inferred from the missing momentum of the event. Validation check: neutral clusters inside a 20-degree cone opposite the K_S0 direction with second moment > 20 cm^2")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on all scan points: real data, inclusive MC, and the matching signal MC
root_files = kskl_alg.execute_on(scan_data + scan_incMC + exMCs)