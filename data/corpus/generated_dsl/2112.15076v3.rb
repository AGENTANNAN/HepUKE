# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation (R-scan continuum: e+e- -> omega pi0 pi0) ###
# R-scan energy points from 2.000 to 3.080 GeV (BOSS 713), given in MeV.
rscan_energies = [2000, 2050, 2100, 2150, 2175, 2200, 2232, 2309,
                  2386, 2396, 2500, 2644, 2646, 2700, 2800,
                  2900, 2950, 2981, 3000, 3020, 3080]

# Real data and inclusive MC for every energy point.
rscan_data  = rscan_energies.map { |e| DatasetManager.real_data.find("713_#{e}") }
rscan_incMC = rscan_energies.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# ConExc "DIY" decay card (mode -2): user-supplied Born cross section (xs_user.txt).
# The virtual photon decays to omega pi0 pi0.  Particle vpho is intentionally
# omitted for a multi-energy scan (injected per energy point by the DSL), and the
# ECMS is likewise handled per point rather than as a single constant.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc -2 vhdr omega pi0 pi0;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k signal-MC events per energy point (same card, one ExclusiveMC per point).
exMCs_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_omega_pi0pi0_rscan"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = '/path/to/xs_user.txt'   # ConExc mode -2 user cross-section file
end

### Event selection (BOSS) ###
alg_name  = "OmegaPi0Pi0"
omega_alg = Algorithm.new(alg_name)
omega_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])

event_selection = Selection.new
  .select_track {                              # exactly one pi+ and one pi-, net charge zero
    cos_theta 0.93                             # |cos(theta)| < 0.93
    Vz        10.0                             # |Vz| < 10 cm
    Vr        1.0                              # Vr < 1 cm
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                             # >= 6 good photons (pi+pi- + six gamma)
    tdc_emc_start     0                        # EMC TDC 0-14
    tdc_emc_end       14
    energyThreshold_b 0.025                    # barrel E > 25 MeV
    energyThreshold_e 0.050                    # endcap E > 50 MeV
    angle_to_track    10.0                     # >= 10 deg from any charged track
    nGam ">=6"
  }
  .pid(method: :probability) {                 # pion PID against kaon and proton
    prob_cut 0.001                             # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]  # pi+ and pi- (charge-conjugation shorthand)
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {    # 1C fit: every gamma-gamma pair -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # best |M(gammagamma) - m_pi0|
    chi2_cut 25
    npi0 ">=3"                                 # one pi0 for omega + two direct pi0
  }
  .kinematic_fit([:pip, :pim, :pi0, :pi0, :pi0]) {   # 4C fit: pi+ pi- pi0 pi0 pi0 hypothesis
    nominal
    constrain_four_momentum
    chi2_cut 200                               # loose; tighter post-fit quality cuts applied later in ROOT
  }

# Best pi0 -> omega (omega -> pi+ pi- pi0) assignment is resolved by the smallest-chi2
# combination iteration inside the nominal 4C fit.
omega_alg
  .note(:omega_pi0_assignment,
        "The pi0 assigned to the omega -> pi+pi-pi0 combination is the one giving M(pi+pi-pi0)
         closest to the omega nominal mass; the nominal 4C fit iterates over the three pi0
         candidates and keeps the combination with the smallest chi2.")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Run over every real-data point, its inclusive MC, and the per-point signal MC.
root_files = omega_alg.execute_on(rscan_data + rscan_incMC + exMCs_signal)