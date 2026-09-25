# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# R-scan energy points spanning 2.000–3.080 GeV (continuum scan).
# Sample names follow the "Rscan_<Ecms(MeV)>" convention for BOTH real data and inclusive MC.
scan_energies = [2000, 2050, 2100, 2125, 2150, 2175, 2200, 2232, 2309, 2386,
                 2396, 2500, 2644, 2646, 2700, 2800, 2900, 2950, 2981, 3000,
                 3020, 3080]

scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("Rscan_#{e}") }   # real-data points
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("Rscan_#{e}") } # matching inclusive MC

# ConExc decay card for e+e- -> K+K-pi0 (continuum / R-scan Born cross-section measurement).
# Mode 8 of the ConExc generator; the pi0 decays via phase space (pi0 -> gamma gamma PHSP).
# `Particle vpho` is deliberately omitted: for a multi-energy scan the DSL injects it
# per energy point automatically.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 K+ K- pi0 ConExc 8;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# A single 100k-event signal sample generated once over all scan points
# (one ExclusiveMC per energy point sharing this same card / cross section).
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "Rscan_KpKmpi0_exclusive_mc"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "RscanKpKmpi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 2.0]})      # scan: the true energy is taken per data point
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                    # charged-track selection
                    cos_theta 0.93                # |cos(theta)| < 0.93
                    Vz        10.0                # |Vz| < 10 cm
                    Vr        1.0                 # Vr < 1 cm in the transverse plane
                    nChrp     ">=1"               # at least one positive track
                    nChrn     ">=1"               # at least one negative track
                    nNet      "==0"               # net charge zero
                  }
                 .select_photon {                 # photon selection
                    energyThreshold_b 0.025       # >= 25 MeV in the EMC barrel
                    energyThreshold_e 0.050       # >= 50 MeV in the EMC endcap
                    nGam   ">=2"                  # at least two photon candidates
                  }
                 .pid(method: :probability) {     # particle identification
                    prob_cut   0.001              # PID probability > 0.001
                    identify :kaon, against: [:pion]  # K+ AND K- (charge-conjugation shorthand), vs pions
                    nkp ">=1"                     # at least one K+
                    nkm ">=1"                     # at least one K-
                  }
                 # Nominal 4C kinematic fit under the K+K-gamma-gamma (K+K-pi0) hypothesis
                 .kinematic_fit([:kp, :km, :gamma, :gamma]) {
                    nominal                       # nominal fit: corrected four-momenta are saved
                    constrain_four_momentum       # 4C energy-momentum constraint
                    chi2_cut 65                   # chi2 < 65
                  }
                 # Competing 4C fit under the K+K-gamma (ISR-phi) hypothesis.
                 # No chi2_cut and no nominal: only the chi2 value is stored so that the
                 # ISR-phi veto (chi2(K+K-pi0) < chi2(K+K-gamma)) can be applied in ROOT.
                 .kinematic_fit([:kp, :km, :gamma]) {
                    constrain_four_momentum
                  }

# Generate the algorithm for the ConExc signal process, then run on all datasets
alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on(scan_data + scan_incMC + exMC_signal)