# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# 15 energy scan points above 4.08 GeV (4.085 - 4.600 GeV); sample name = [BOSS]_[ECMS(MeV)]
sample_names = %w[
  703_4090
  703_4180
  703_4190
  703_4200
  703_4210
  703_4220
  703_4230
  703_4237
  703_4246
  703_4260
  703_4270
  703_4280
  703_4360
  703_4420
  703_4600
]
data_points = sample_names.map { |n| DatasetManager.real_data.find(n) }      # real data at each point
incMC_points = sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # matching inclusive MC

# Decay card for e+e- -> pi+ pi- D0 Dbar0, with D0 -> K- pi+ (charge conjugate implied)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC (200k events) of the pi+ pi- D0 Dbar0 (D0 -> K- pi+) phase-space mode,
# generated separately for every energy point.
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipiD0D0bar"   # auto-suffixed per energy point
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pipiD0D0bar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})   # beam energy (per-point value handled by the job)
            .set_alias({"std::vector<double>" => "Vdouble"})

# One common selection chain shared by all energy points
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta   0.93      # |cos(theta)| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # Vr < 1 cm
                  nChrp       ">=2"     # at least 2 positive tracks
                  nChrn       ">=2"     # at least 2 negative tracks
                }
               .select_photon {         # Photon selection (shower quality)
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    20.0   # angle to nearest charged track > 20 deg
                  energyThreshold_b 0.025  # 25 MeV in the barrel
                  energyThreshold_e 0.050  # 50 MeV in the endcap
                }
               .pid(method: :probability) {  # PID: at least one K+ and one K-
                  prob_cut 0.001
                  identify :kaon, against: [:pion]   # K+ and K- vs pions
                  nkp ">=1"
                  nkm ">=1"
                }
               .remove([:kp <= :chrgp, :km <= :chrgn])     # take identified kaons out
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks are pi+ / pi-
               .kinematic_fit([:pip, :pip, :pim, :pim, :kp, :km]) {  # 4C fit on the pi pi K K system
                  nominal
                  constrain_four_momentum
                  chi2_cut 56
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs)