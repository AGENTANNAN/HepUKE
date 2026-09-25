# =============================================================================
# e+e- -> eta pi+ pi- (eta -> gamma gamma) at the 713 R-scan energy points
# (continuum measurement -> ConExc generator, ISR included)
# =============================================================================

### Dataset preparation ###
# 713 R-scan points between 2.000 and 3.080 GeV (sample name convention <BOSS>_<sample>)
rscan_points = %w[2000 2050 2100 2150 2175 2200 2232 2309 2386 2396 2500 2644 2646 2700 2800 2900 2950 2981 3000 3020 3080]

rscan_data  = rscan_points.map { |e| DatasetManager.real_data.find("713_Rscan_#{e}") }     # real data at each point
rscan_incMC = rscan_points.map { |e| DatasetManager.inclusive_mc.find("713_Rscan_#{e}") } # matching inclusive MC at each point

# ConExc decay card (continuum / R-scan): mode 37 = e+e- -> eta pi+ pi-, eta -> gamma gamma, ISR included.
# No KKMC and no "Particle vpho" line: the DSL injects Particle vpho <ECMS> 0.0 per energy point.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 37;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC, one sample per R-scan energy point (shared card, ISR modelled by ConExc)
exMCs = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_rscan_eta_pipi"  # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtaPiPi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.080]}) # representative CMS energy of the scan

event_selection = Selection.new
  .select_track {              # exactly two charged tracks, net charge 0
    cos_theta 0.93             # |cos(theta)| < 0.93 (polar angle of charged tracks)
    Vz        10.0             # |Vz| < 10 cm (vertex along beam axis)
    Vr        1.0              # Vr < 1 cm (vertex in transverse plane)
    nChrp     "==1"            # exactly one positive track
    nChrn     "==1"            # exactly one negative track
    nNet      "==0"            # net charge zero
  }
  .select_photon {             # at least two good photons
    tdc_emc_start     0        # TDC-EMC window start
    tdc_emc_end       14       # TDC-EMC window end
    energyThreshold_b 0.025    # 25 MeV minimum shower energy in the barrel
    energyThreshold_e 0.050    # 50 MeV minimum shower energy in the endcap
    angle_to_track    10.0     # photon-track opening angle > 10 degrees
    nGam              ">=2"    # at least two photons (eta -> gamma gamma)
  }
  .pid(method: :probability) { # probability-method PID: both tracks are pions
    prob_cut 0.001             # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K / p hypotheses
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {  # 4C fit: e+e- -> pi+ pi- gamma gamma
    nominal                    # this is the nominal fit (corrected four-momenta are kept)
    constrain_four_momentum    # constrain the total four-momentum to the CMS system
    chi2_cut 100               # chi^2 < 100
  }

# BOSS-side procedures that have no formal DSL construct
algorithm
  .note(:background_veto, "Bhabha (e+e- -> e+e-) events suppressed by requiring E/p < 0.8 for each charged track before the 4C kinematic fit")
  .note(:helicity_angle_cut, "each photon required to satisfy cos(theta_gamma) < 0.95 in the eta helicity frame")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC (one job per energy point)
root_files = algorithm.execute_on(rscan_data + rscan_incMC + exMCs)