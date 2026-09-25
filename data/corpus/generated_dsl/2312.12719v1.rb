# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# 12 ISR scan energy points taken with BOSS 703 ...
isr_energies_703 = %w[4009 4180 4190 4200 4210 4220 4230 4260 4270 4360 4420 4600]
isr_data  = isr_energies_703.map { |e| DatasetManager.real_data.find("703_#{e}") }      # real data at each ISR scan point
isr_incMC = isr_energies_703.map { |e| DatasetManager.inclusive_mc.find("703_#{e}") }  # inclusive MC at each ISR scan point

# ... plus the dedicated 3.773 GeV point taken with BOSS 712
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

all_data  = isr_data  + [data_3773]
all_incMC = isr_incMC + [incMC_3773]

# Decay card for the signal: e+e- -> (untagged ISR) Sigma+ Sigma-,
# Sigma+ -> p+ pi0, anti-Sigma- -> anti-p- pi0, pi0 -> gamma gamma (all phase space)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# One exclusive signal MC sample per scan point, 100k events each
exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_sigmapm_isr"   # auto-suffixed per scan point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "SigmapmISR"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})           # nominal CMS energy of the psi(4260) sample
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                     # exactly two good charged tracks, net charge zero
      cos_theta   0.93                # |cos(theta)| < 0.93
      Vz          10.0                # |Vz| < 10 cm
      Vr          1.0                 # Vr < 1 cm
      nChrp       "==1"               # one positive track (=> at least one + and one -)
      nChrn       "==1"               # one negative track
      nNet        "==0"               # net charge zero
  }
  .select_photon {                    # at least four photons (two pi0 candidates)
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0          # >= 10 deg from any charged track
      energyThreshold_b 0.025         # 25 MeV in the EMC barrel
      energyThreshold_e 0.050         # 50 MeV in the EMC endcap
      nGam              ">=4"         # at least 4 photons
  }
  .pid(method: :probability) {        # proton / anti-proton identification
      prob_cut 0.001
      identify :proton, against: [:pion, :kaon]   # p+ and anti-p- vs pi+ and K+
      nprp ">=1"
      nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # 1C fit to reconstruct pi0 from gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=2"                                  # two pi0 candidates (one per Sigma)
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {      # 4C fit on the p pbar pi0 pi0 system
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

# Attach the decay card and render the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Run on real data, inclusive MC and exclusive signal MC for every scan point
root_files = my_algorithm.execute_on(all_data + all_incMC + exMC_signal)