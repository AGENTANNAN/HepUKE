# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# BOSS 703 XYZ scan energy points spanning 3.808-4.600 GeV
data_energy_points = %w[
  703_3810 703_3872 703_3900 703_4009 703_4090 703_4180 703_4190 703_4200
  703_4210 703_4220 703_4230 703_4237 703_4245 703_4246 703_4260 703_4270
  703_4280 703_4310 703_4360 703_4390 703_4420 703_4470 703_4530 703_4575 703_4600
]
data_points = data_energy_points.map { |name| DatasetManager.real_data.find(name) }

# Matched inclusive MC (available for the subset of scan points in this range)
incMC_energy_points = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230 703_4237
  703_4246 703_4260 703_4270 703_4280 703_4360 703_4420 703_4600
]
incMC_points = incMC_energy_points.map { |name| DatasetManager.inclusive_mc.find(name) }

### Decay cards ###
# ConExc is used for the ISR treatment of the continuum production of pi0 pi0 J/psi.
# The DSL auto-detects the literal token ConExc and injects "Particle vpho <ECMS> 0.0"
# per energy point, so Particle vpho is deliberately omitted for this multi-energy scan.
# J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay vpho
    1.000  ConExc  74110;
    Enddecay

    Decay J/psi
    1.000  e+  e-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay vpho
    1.000  ConExc  74110;
    Enddecay

    Decay J/psi
    1.000  mu+  mu-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC (100k events per scan point, one sample per J/psi leptonic mode) ###
exMCs_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pi0pi0Jpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMCs_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pi0pi0Jpsi_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0Pi0Jpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 4.26] })
            .set_alias({ "std::vector<double>" => "Vdouble" })

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93      # |cos(theta)| < 0.93
                  Vz        10.0      # |Vz| < 10 cm
                  Vr        1.0       # Vr < 1 cm
                  nChrp     "==1"     # exactly one positive track
                  nChrn     "==1"     # exactly one negative track
                  nNet      "==0"     # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0    # > 10 degrees from any charged track
                  energyThreshold_b 0.025   # 25 MeV in the barrel
                  energyThreshold_e 0.050   # 50 MeV in the endcap
                  nGam              ">=4"   # pi0 pi0 -> 4 photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  # E/p-based lepton ID (e: E/p>0.7, mu: E/p<0.3) is not expressible;
                  # the high-momentum lepton helper is the closest DSL approximation,
                  # giving the combined l+/l- lists used in the fit below.
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon, :proton]  # separate pions from K and p
                  nlp  "==1"   # one l+
                  nlm  "==1"   # one l-
                }
               # Reconstruct pi0 candidates from gamma gamma with mass-constrained Kalman fits (chi2 < 25)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=2"   # at least two pi0 candidates (one per pi0 of the final state)
                }
               # 4C kinematic fit to l+ l- pi0 pi0; the best pi0 pi0 pairing (out of the
               # three possible pairings) is chosen automatically by smallest chi2
               .kinematic_fit([:lp, :lm, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 75
                }

my_algorithm
  .note(:pid_correction_method, "lepton identification assigns the positive/negative tracks as l+/l- using E/p (electron E/p > 0.7, muon E/p < 0.3); this E/p criterion has no DSL method, so identify_high_momentum_leptons (track momentum / EMC energy based) is used as the closest available approximation")
  .note(:isr_and_vacuum_polarization, "ISR up to second order and the measured sigma0(m) are supplied by the ConExc generator via the Decay vpho block; the vacuum-polarization correction factors are taken from the ConExc generator log for the Born cross-section extraction")
  .with_decay_card(decay_card_ee)
  .apply(event_selection)

# Execute on all scan points, their matched inclusive MC, and both signal MC sets
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_ee + exMCs_mumu)