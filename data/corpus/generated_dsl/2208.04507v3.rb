# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# e+e- -> omega pi+ pi- measured at c.m. energies 2.000-3.080 GeV over 19 points
# (713 R-scan, 647 pb^-1 total). Real data and inclusive MC at the same points.
scan_energies = [2000, 2050, 2100, 2150, 2175, 2200, 2232, 2309, 2386, 2396,
                 2500, 2644, 2646, 2900, 2950, 2981, 3000, 3020, 3080]
scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("713_#{e}") }
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# Decay card (EvtGen): e+e- -> omega pi+pi- (phase space),
# omega -> pi+ pi- pi0 (Dalitz), pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0  omega  pi+  pi-    PHSP;
    Enddecay

    Decay omega
    1.0  pi+  pi-  pi0    OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0  gamma  gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: the same e+e- -> omega pi+pi- sample is generated at every
# scan energy point (one related dataset per point), 500k events each.
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_omega_pi_pi"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaPiPiPi"
omega_algorithm = Algorithm.new(alg_name)
omega_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
               .set_constant({"ECMS" => [:double, 3.08]})
               .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    # Four charged tracks, net charge zero; each track |cos(theta)| < 0.93,
    # |Vz| < 10 cm, Vr < 1 cm
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
    }
    # At least two photons: E > 25 MeV (barrel) / > 50 MeV (endcap),
    # EMC time within 0-14 (0-700 ns), and >= 10 degrees from any charged track
    .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=2"
    }
    # Probability-based PID (prob > 0.001): identify pions against kaons;
    # require two pi+ and two pi- (pi+ pi- pi+ pi-)
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        npip "==2"
        npim "==2"
    }
    .assign({:chrgp => :pip, :chrgn => :pim})
    # 4C kinematic fit to e+e- -> 2(pi+pi-) gamma gamma with four-momentum
    # conservation; pi0 built from the gamma gamma pair (invariant mass 0.120-0.150 GeV/c^2)
    .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) {
        nominal
        invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)
        constrain_four_momentum
        chi2_cut 200
    }

omega_algorithm
    .note(:beam_energy_scan, "signal is measured over 19 c.m. energy points (2.000-3.080 GeV, 647 pb^-1 total); the single ECMS constant is a placeholder - the true beam energy differs per scan point and is set through the related dataset of the exclusive MC, so the CMS four-momentum in the kinematic fit must be built from the per-run measured beam energy")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on all 19 scan points: real data + inclusive MC + signal exclusive MC
root_files = omega_algorithm.execute_on(scan_data + scan_incMC + exMCs_signal)