# =============================================================================
#  e+e- -> omega pi0   and   e+e- -> omega eta
#  34 energy points in the range 3.773 - 4.701 GeV
#  omega -> pi+ pi- pi0 ,  pi0 -> gamma gamma ,  eta -> gamma gamma
# =============================================================================

### Dataset preparation ###

# The 34 real-data energy points spanning 3.773 - 4.701 GeV
data_points = %w[
  712_3773 712_3780 712_3800 712_3815 712_3825
  712_3835 712_3845 712_3860 712_3885
  703_3900 703_4009 703_4090 705_4130 705_4160
  703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4246 703_4260 703_4270
  703_4280 705_4290 703_4310 705_4340 703_4360
  703_4390 703_4420 706_4610 706_4660 706_4700
].map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC is only available for a subset of the scan points
incMC_points = %w[
  705_4130 705_4160
  703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4246 703_4260 703_4270
  703_4280 703_4360 703_4420
  706_4610 706_4660 706_4700
].map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay card for e+e- -> omega pi0 (continuum process, KKMC top mother convention)
decay_card_omega_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 omega pi0 PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for e+e- -> omega eta
decay_card_omega_eta = <<~DECAYCARD
  Decay psi(4260)
  1.0000 omega eta PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for each channel, one sample per energy point of the scan
exMCs_omega_pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_pi0"
  config.events        = 100000
  config.decay_card    = decay_card_omega_pi0
  config.cross_section = :default
end

exMCs_omega_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_eta"
  config.events        = 100000
  config.decay_card    = decay_card_omega_eta
  config.cross_section = :default
end

### Event selection (BOSS) ###

# Common selection chain shared by the two channels (they differ only in the
# final-state particle assignment of the kinematic fit).
event_selection_common = Selection.new
  .select_track {                       # exactly one pi+ and one pi-
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"                     # one positively charged track
    nChrn     "==1"                     # one negatively charged track
    nTot      "==2"                     # two charged tracks in total
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # at least four good photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0              # at least 10 degrees from any charged track
    energyThreshold_b 0.025             # barrel  > 25 MeV
    energyThreshold_e 0.050             # endcap  > 50 MeV
    nGam              ">=4"
  }
  .pid(method: :probability) {          # probability method PID
    prob_cut 0.001                      # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]   # pi+ / pi- separated from K and p
    npip "==1"
    npim "==1"
  }

# -----------------------------------------------------------------------------
# Channel I : e+e- -> omega pi0  ( omega -> pi+ pi- pi0, both pi0 -> gamma gamma )
# -----------------------------------------------------------------------------
omega_pi0_selection = event_selection_common.dup
  # mass-constrained Kalman fits on gamma-gamma pairs:
  # the pi0 from the omega and the bachelor pi0
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # 6C fit: 4-momentum conservation + the two nominal pi0 mass constraints
  .kinematic_fit([:pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # competing 5C fit: bachelor pi0 left without the nominal mass constraint
  # (it enters as its two free photons, while the omega-side pi0 stays constrained)
  .kinematic_fit([:pip, :pim, :pi0, :gamma, :gamma]) {
    constrain_four_momentum
    chi2_cut 200
  }

alg_name_omega_pi0 = "OmegaPi0"
alg_omega_pi0 = Algorithm.new(alg_name_omega_pi0)
alg_omega_pi0.set_header(["#{alg_name_omega_pi0}Alg/#{alg_name_omega_pi0}.h"])
             .set_constant({"ECMS" => [:double, 4.178]})   # representative scan energy
             .set_alias({"std::vector<double>" => "Vdouble"})
# ECMS is a single constant while the analysis runs over 34 beam energies; the
# per-run beam energy is taken from the conditions DB at job time.
alg_omega_pi0.note(:beam_energy_per_point,
  "dataset spans 34 energy points from 3.773 to 4.701 GeV; the CMS four-momentum " \
  "is built from the per-run beam energy read from the conditions DB, the ECMS " \
  "constant is only a representative fallback")

alg_omega_pi0.with_decay_card(decay_card_omega_pi0).apply(omega_pi0_selection)
root_files_omega_pi0 = alg_omega_pi0.execute_on(data_points + incMC_points + exMCs_omega_pi0)

# -----------------------------------------------------------------------------
# Channel II : e+e- -> omega eta  ( omega -> pi+ pi- pi0, pi0/eta -> gamma gamma )
# -----------------------------------------------------------------------------
omega_eta_selection = event_selection_common.dup
  # reconstruct the pi0 from the omega
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # reconstruct the eta
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # 6C fit: 4-momentum conservation + the nominal pi0 and eta mass constraints
  .kinematic_fit([:pip, :pim, :pi0, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # competing 5C fit: eta left without the nominal mass constraint
  # (enters as its two free photons, the pi0 from the omega stays constrained)
  .kinematic_fit([:pip, :pim, :pi0, :gamma, :gamma]) {
    constrain_four_momentum
    chi2_cut 200
  }

alg_name_omega_eta = "OmegaEta"
alg_omega_eta = Algorithm.new(alg_name_omega_eta)
alg_omega_eta.set_header(["#{alg_name_omega_eta}Alg/#{alg_name_omega_eta}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})   # representative scan energy
            .set_alias({"std::vector<double>" => "Vdouble"})
alg_omega_eta.note(:beam_energy_per_point,
  "dataset spans 34 energy points from 3.773 to 4.701 GeV; the CMS four-momentum " \
  "is built from the per-run beam energy read from the conditions DB, the ECMS " \
  "constant is only a representative fallback")

alg_omega_eta.with_decay_card(decay_card_omega_eta).apply(omega_eta_selection)
root_files_omega_eta = alg_omega_eta.execute_on(data_points + incMC_points + exMCs_omega_eta)