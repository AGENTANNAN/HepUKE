# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Representative scan points of the e+e- -> Sigma- anti-Sigma+ measurement (sqrt(s) = 3.51-4.95 GeV).
# Sample names follow the BESIII convention [BOSS version]_[CMS energy in MeV].
data_3773 = DatasetManager.real_data.find("712_3773")   # psi(3770) point
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

scan_data = [data_3773, data_4600, data_4610, data_4620, data_4640, data_4660,
             data_4680, data_4700, data_4740, data_4750, data_4780, data_4840,
             data_4914, data_4946]

# Matching inclusive MC samples at the same scan points.
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

scan_incMC = [incMC_3773, incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660,
              incMC_4680, incMC_4700, incMC_4740, incMC_4750, incMC_4780, incMC_4840,
              incMC_4914, incMC_4946]

# Decay card for the signal process: e+e- -> Sigma- anti-Sigma+, Sigma- -> n pi-,
# anti-Sigma+ -> anti-n pi+  (EvtGen syntax; psi(4260) is the KKMC top-mother convention).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Sigma- anti-Sigma+ PHSP;
  Enddecay

  Decay Sigma-
  1.000 n0 pi- PHSP;
  Enddecay

  Decay anti-Sigma+
  1.000 anti-n0 pi+ PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC at every scan point (same decay card, cross section, statistics).
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_sig_SigmamSigmabar"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "SigmamSigmabar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # nominal energy of the analysis

event_selection = Selection.new
event_selection
  .select_track {                 # charged-track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        30.0                # |Vz| < 30 cm
    Vr        10.0                # |Vr| < 10 cm
    nChrp     "==1"               # exactly one positively charged track
    nChrn     "==1"               # exactly one negatively charged track
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # EMC shower selection
    tdc_emc_start     0           # EMC time window [0, 700] ns
    tdc_emc_end       14
    angle_to_track    15.0        # > 15 deg w.r.t. any charged track
    energyThreshold_b 0.025       # barrel: E > 25 MeV
    energyThreshold_e 0.050       # endcap: E > 50 MeV
    nGam              ">=1"       # at least one shower
  }
  .pid(method: :probability) {    # PID: pions separated from kaons/protons
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- via charge-conjugation shorthand
    npip "==1"
    npim "==1"
  }
  # The antineutron is reconstructed from the highest-energy EMC shower; it is carried by the
  # EMC-shower (photon) candidate list (shower-shape criteria are captured in the notes).
  .for_each(:gamma) do
    where { energy > 0.62 }       # energy-dependent threshold: 0.62 GeV at 3.773 GeV
    best { maximize { energy } }  # pick the highest-energy EMC shower
    remove                        # keep only that antineutron candidate
  end
  # 1-C kinematic fit: constrain the (pi+ nbar) invariant mass to the nominal anti-Sigma+ mass.
  .kinematic_fit([:pip, :pim, :gamma]) do
    nominal
    invariant_mass_of(:pip, :gamma).constrain_to_nominal_mass_of(:"anti-Sigma+")
    chi2_cut 60
  end

my_algorithm
  .note(:nbar_reconstruction, "the antineutron is the highest-energy EMC shower with deposited
    energy above an energy-dependent threshold (0.62 GeV at 3.773 GeV), second moment S > 18 cm^2,
    more than 16 hits and an angle > 15 deg to any charged track; the threshold is re-evaluated per
    scan point")
  .note(:recoil_mass_window, "partial-reconstruction/recoil technique: the neutron is left
    undetected and the Sigma- is inferred from the recoil against the reconstructed
    anti-Sigma+ (pi+ nbar); events are required to have the recoil mass M(pi+ nbar) in
    [1.10, 1.35] GeV/c^2 (Sigma- mass region)")
  .note(:nbar_direction_constraint, "the 1-C kinematic fit additionally constrains the
    antineutron direction to the EMC shower position")
  .note(:tof_veto, "cosmic-ray events are rejected by requiring the MDC time-of-flight
    difference between the pi+ and the pi- to be greater than -5 ns")
  .note(:isr_correction, "detection efficiencies and ISR corrections are obtained iteratively
    from the exclusive signal MC; the ISR/continuum treatment is applied when extracting the
    Born cross sections and in the psi(3770) -> Sigma- anti-Sigma+ search")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on all real data, inclusive MC and signal MC samples.
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMCs_signal)