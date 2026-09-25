# frozen_string_literal: true
# BOSS event-selection DSL for the measurement of e+e- -> K_S0 K_L0 pi0 cross sections
# over the 19 BOSS 713 R-scan energy points spanning 2.000 - 3.080 GeV.

### Dataset preparation ###
# BOSS 713 R-scan energy points (2.000 - 3.080 GeV)
rscan_energies = %w[
  2000 2050 2100 2150 2175 2200 2232 2309 2386 2396
  2500 2700 2800 2900 2950 2981 3000 3020 3080
]
# Real data and inclusive MC at each energy point (sample name = "<BOSS>_<SampleName>")
rscan_data  = rscan_energies.map { |e| DatasetManager.real_data.find("713_Rscan_#{e}") }
rscan_incMC = rscan_energies.map { |e| DatasetManager.inclusive_mc.find("713_Rscan_#{e}") }

# ConExc decay card: continuum / R-scan generator modelling ISR up to second order and
# vacuum polarisation. e+e- -> K_S0 K_S0 pi0 (ConExc mode 9); one K_S0 is treated as the
# missing K_L0 in the analysis. No "Particle vpho" line -- the DSL injects it per energy
# point automatically for multi-energy scans.
conexc_decay_card = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 9;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC at every energy point (same card, one MC per point)
exMC_conexc = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_conexc_ksklpi0"   # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = conexc_decay_card
  config.cross_section = :default
end
exMC_conexc.each { |m| m.save_to_config(format: :yaml, file_path: 'conexc_mc_config') }

### Event selection (BOSS) ###
alg_name = "KsKlPi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.080]})   # per-run value in the R scan
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  # Charged-track selection: exactly one pi+ and one pi-
  .select_track {
    cos_theta 0.93       # |cos(theta)| < 0.93
    Vz        10.0       # |Vz| < 10 cm
    Vr        1.0        # Vr < 1 cm
    nChrp     "==1"      # exactly one positive track
    nChrn     "==1"      # exactly one negative track
    nNet      "==0"      # net charge zero
  }
  # Photon selection
  .select_photon {
    tdc_emc_start     0      # EMC TDC start = 0
    tdc_emc_end       14     # EMC TDC end   = 14
    energyThreshold_b 0.025  # 25 MeV barrel threshold
    energyThreshold_e 0.050  # 50 MeV endcap threshold
    angle_to_track    10.0   # > 10 deg to nearest charged track
    nGam              ">=2"  # at least two photons
  }
  # PID: pi+ / pi-
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  # K_S0 <- pi+ pi- through a secondary-vertex fit; pick the combination closest to the
  # nominal K_S0 mass and drop the used tracks
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # pi0 <- gamma gamma (Kalman 1C mass-constrained fit), chi2 < 30, at least one pi0
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    npi0 ">=1"
  }
  # ISR-veto competing-hypothesis fit: e+e- -> gamma K_S0 K_L0 (chi2 stored for ROOT veto)
  .kinematic_fit([:gamma, :K_S0, :K_L0]) {
    miss_track_of :K_L0
    constrain_four_momentum
  }
  # Nominal 3C fit: K_L0 missing, K_S0 mass constrained to nominal, 4-momentum conservation
  .kinematic_fit([:K_S0, :pi0, :K_L0]) {
    nominal
    miss_track_of :K_L0
    constrain_four_momentum
    chi2_cut 200
  }

my_algorithm
  .note(:background_veto, "ISR veto: a 1C fit is performed under the e+e- -> gamma K_S0 K_L0 hypothesis alongside the nominal 3C fit; the competing chi2 is stored unconditionally and the veto (chi2_3C(K_S0 K_L0 pi0) < chi2_1C(gamma K_S0 K_L0)) is applied at the ROOT level.")
  .note(:efficiency_curve, "K_S0 and pi0 candidate windows are applied to the measured (pre-nominal-fit) quantities: |M(pi+pi-) - m_K_S0| < 12 MeV, |M(gamma gamma) - m_pi0| < 0.015 GeV, and K_S0 decay length > 2 sigma from the secondary-vertex fit; these enter the efficiency correction.")
  .note(:signal_mc_treatment, "signal MC is generated as e+e- -> K_S0 K_S0 pi0 (ConExc mode 9) and one K_S0 is treated as the missing K_L0 in the cross-section extraction; ISR and vacuum-polarisation corrections are taken from the ConExc generator output.")
  .with_decay_card(conexc_decay_card)
  .apply(event_selection)

# Execute on real data, inclusive MC and the per-energy-point exclusive signal MC
root_files = my_algorithm.execute_on(rscan_data + rscan_incMC + exMC_conexc)