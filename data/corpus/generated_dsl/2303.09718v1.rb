# ==============================================================================
# Continuum e+e- -> pi+ pi- omega : omega -> pi+ pi- pi0 (Dalitz), pi0 -> gamma gamma
# at 24 CMS energies (4.009 - 4.680 GeV), for a helicity-amplitude analysis of
# f0(500), f0(980), f2(1270), f0(1370), b1(1235)+- and rho(1450)+-
# fitted simultaneously over the energy points.
# Continuum (ISR) production is modelled with the ConExc generator.
# Charge conjugation is implied throughout.
# ==============================================================================

### Dataset preparation ###
# 24 real-data samples covering sqrt(s) = 4.009 - 4.680 GeV
# (BOSS version per the BES3 dataset table; 703 below 4.60 GeV, 706 above)
data_points = [
  DatasetManager.real_data.find("703_4009"),   # 4.009 GeV
  DatasetManager.real_data.find("703_4180"),   # 4.180 GeV
  DatasetManager.real_data.find("703_4190"),   # 4.190 GeV
  DatasetManager.real_data.find("703_4200"),   # 4.200 GeV
  DatasetManager.real_data.find("703_4210"),   # 4.210 GeV
  DatasetManager.real_data.find("703_4220"),   # 4.220 GeV
  DatasetManager.real_data.find("703_4230"),   # 4.230 GeV
  DatasetManager.real_data.find("703_4237"),   # 4.237 GeV
  DatasetManager.real_data.find("703_4246"),   # 4.246 GeV
  DatasetManager.real_data.find("703_4260"),   # 4.260 GeV
  DatasetManager.real_data.find("703_4270"),   # 4.270 GeV
  DatasetManager.real_data.find("703_4280"),   # 4.280 GeV
  DatasetManager.real_data.find("703_4310"),   # 4.310 GeV
  DatasetManager.real_data.find("703_4360"),   # 4.360 GeV
  DatasetManager.real_data.find("703_4390"),   # 4.390 GeV
  DatasetManager.real_data.find("703_4420"),   # 4.420 GeV
  DatasetManager.real_data.find("703_4470"),   # 4.470 GeV
  DatasetManager.real_data.find("703_4530"),   # 4.530 GeV
  DatasetManager.real_data.find("703_4575"),   # 4.575 GeV
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV
  DatasetManager.real_data.find("706_4610"),   # 4.610 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.620 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.640 GeV
  DatasetManager.real_data.find("706_4680")    # 4.680 GeV
]

# Inclusive MC (background) only at 4.180 and 4.260 GeV
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# Decay card for the continuum signal process, using the ConExc generator
# (ISR + measured Born cross section). The literal token "ConExc" switches the
# DSL to the no-KKMC template and auto-injects `Particle vpho <ECMS>` per energy
# point, so no `Particle vpho` line is written for this multi-energy scan.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 pi+ pi- omega ConExc 4;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive ConExc MC per energy point (same card / cross section)
exMC_signals = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_conexc_omegapippi"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaPiPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                         # Charged-track selection
      cos_theta 0.93                      # |cos(theta)| < 0.93
      Vz        10.0                      # |Vz| < 10 cm
      Vr        1.0                       # Vr < 1 cm
      nChrp     ">=2"                     # at least two positive tracks
      nChrn     ">=2"                     # at least two negative tracks
      nNet      "==0"                     # net charge zero
  }
  .select_photon {                        # Photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0              # > 10 deg from any charged track
      energyThreshold_b 0.025             # > 25 MeV in the barrel
      energyThreshold_e 0.050             # > 50 MeV in the endcap
      nGam              ">=2"             # at least two photons
  }
  .pid(method: :probability) {            # Pion identification
      prob_cut 0.001                      # PID probability > 0.001
      identify :pion, against: [:kaon, :proton]  # pi+/pi- vs K and p (both charges)
      npip ">=2"                          # at least two pi+
      npim ">=2"                          # at least two pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from photon pairs (1C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"                          # at least one pi0 candidate
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0]) {   # 5C nominal fit (4C + pi0 mass)
      nominal
      constrain_four_momentum
      chi2_cut 60
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :gamma]) {  # competing pi+pi-pi+pi-pi0 gamma hypothesis
      constrain_four_momentum             # store chi2 for ROOT-level background suppression
  }

# Post-5C-fit selection criteria that have no BOSS-level DSL expression are kept as notes.
my_algorithm
  .note(:e_over_p_cut, "E_EMC/p < 0.9 applied to pion candidates not originating from omega decays (suppresses non-omega multi-pion background); applied after the 5C fit")
  .note(:background_veto, "veto |M(pi+pi-) - M(K_S0)| in (0.49, 0.51) GeV/c^2 (K_S0 -> pi+pi-) and |M(pi+pi-pi0) - M(chi_c0)| in (3.39, 3.44) GeV/c^2 (chi_c0 background), applied after the 5C fit")
  .note(:omega_candidate_selection, "when several omega candidates pass, the one whose M(pi+pi-pi0) is closest to the nominal omega mass is selected")
  .note(:signal_region, "omega signal region (0.76, 0.82) GeV/c^2; sidebands (0.68, 0.74) and (0.84, 0.90) GeV/c^2 used for background estimation in the simultaneous amplitude fit")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on the 24 data points, the two inclusive MC samples and the per-energy ConExc MC
root_files = my_algorithm.execute_on(data_points + [incMC_4180, incMC_4260] + exMC_signals)