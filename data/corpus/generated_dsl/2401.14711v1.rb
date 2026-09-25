# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# R-scan energy points (2015 R-scan restored in BOSS 713) covering 2.00 – 3.08 GeV
scan_energies = [3080, 3020, 3000, 2981, 2950, 2900,          # high-mass points
                 2644, 2646, 2500, 2396, 2386, 2309, 2232,    # intermediate points
                 2200, 2175, 2150, 2100, 2050, 2000]          # 19 CM energies in total

# Real data and corresponding inclusive MC samples, one per CM energy point
scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("713_#{e}") }
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# ConExc decay card for the continuum process e+e- -> pi+ pi- pi0 (mode 7).
# The literal token "ConExc" makes the DSL use the no-KKMC simulation template and
# inject "Particle vpho <ECMS> 0.0" for each energy point, so no Particle vpho line
# is written here.  ConExc provides the Born cross section and ISR up to 2nd order,
# and its generator log supplies the ISR correction factor.
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc 7;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 50k events per scan point, generated with ConExc
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
    config.sample_name   = "exmc_conexc_pipimpimpi0"  # auto-suffixed per energy point
    config.events        = 50_000
    config.decay_card    = decay_card_conexc
    config.cross_section = :default                  # ConExc Born cross section model
end

### Event selection (BOSS) ###
alg_name     = "PipPimPi0Rscan"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.08]})  # reference value; supplied per scan point at job time
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                       # Charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        1.0                     # |Vr| < 1 cm
        nChrp     ">=1"                   # at least one positive track
        nChrn     ">=1"                   # at least one negative track
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # Photon selection
        tdc_emc_start     0               # EMC TDC start
        tdc_emc_end       14              # EMC TDC end
        angle_to_track    10.0            # > 10 deg to nearest charged track
        energyThreshold_b 0.025           # barrel energy > 25 MeV
        energyThreshold_e 0.050           # endcap energy > 50 MeV
        nGam              ">=2"           # at least two photons (for pi0 -> gamma gamma)
    }
    .pid(method: :probability) {          # PID: probability method
        prob_cut 0.001                    # PID probability > 0.001
        identify :pion, against: [:kaon, :proton]  # pi+ and pi- identified together
        npip     ">=1"                    # at least one pi+
        npim     ">=1"                    # at least one pi-
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1-C Kalman fit: pi0 -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                       # chi2 < 25
        npi0     ">=1"                    # at least one pi0 candidate
    }
    .kinematic_fit([:pip, :pim, :pi0]) {  # Final 4C kinematic fit on pi+ pi- pi0
        nominal                           # nominal fit: corrected four-momenta are saved
        constrain_four_momentum           # constrain total four-momentum to the CMS energy
        chi2_cut 50                       # chi2 < 50
    }

# Generate the algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_conexc).apply(event_selection)

# Execute on all scan points: real data, inclusive MC and ConExc signal MC
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMCs_signal)