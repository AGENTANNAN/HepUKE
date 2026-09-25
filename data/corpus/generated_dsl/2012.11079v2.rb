# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Study of e+e- -> p pbar p pbar over the c.m. energy scan 4.009-4.600 GeV.
# BOSS 703 = XYZ scan; sample names follow [BOSS version]_[CMS energy in MeV].
real_data_names = %w[
  703_4009 703_4090 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4237 703_4245 703_4246 703_4260
  703_4270 703_4280 703_4310 703_4360 703_4390 703_4420
  703_4470 703_4530 703_4575 703_4600
]
data_points = real_data_names.map { |n| DatasetManager.real_data.find(n) }

# Matching inclusive MC samples at the scan points where they are available.
incMC_names = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4246 703_4260 703_4270 703_4280
  703_4360 703_4420 703_4600
]
incMC_points = incMC_names.map { |n| DatasetManager.inclusive_mc.find(n) }

# Decay card for the signal e+e- -> 2(p pbar): psi(4260) (KKMC top-mother
# convention) decaying to two proton-antiproton pairs in phase space.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 p+ anti-p- p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC generated at every energy point of the scan
# (same decay card / cross section, only the related dataset differs).
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ppbar_ppbar"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PPbarPPbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})       # reference c.m. energy (scan span 4.009-4.600 GeV)
            .set_alias({"std::vector<double>" => "Vdouble"})

# Single selection chain shared by real data, inclusive MC and exclusive signal MC.
event_selection = Selection.new
event_selection
  .select_track {
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr        1.0      # Vr < 1 cm
      nChrp     "==2"    # exactly two positive charged tracks
      nChrn     "==2"    # exactly two negative charged tracks
      nNet      "==0"    # net charge zero
      # no photon requirement in this analysis
  }
  .pid(method: :probability) {
      prob_cut 0.001                              # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]   # p+ and anti-p- separated from K and pi
      nprp "==2"                                  # exactly two protons
      nprm "==2"                                  # exactly two antiprotons
  }
  # 3C kinematic fit over the four (anti)proton tracks, constraining the
  # total three-momentum to the c.m. system.
  .kinematic_fit([:prp, :prm, :prp, :prm]) {
      nominal                    # nominal fit
      constrain_three_momentum   # 3C constraint
      chi2_cut 60                # chi^2 < 60
  }

my_algorithm
  .note(:energy_scan_ecms, "the analysis spans 22/23 c.m. energy points between 4.009 and 4.600 GeV; a single ECMS constant cannot encode the per-point beam energy, so each dataset must be processed with its own measured c.m. energy")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the same selection on real data, inclusive MC and exclusive signal MC.
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_signal)