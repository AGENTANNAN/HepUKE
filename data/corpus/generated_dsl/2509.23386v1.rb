### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Signal decay card: ψ(3686) → γ χ_cJ (J=0,1,2) with equal weight, χ_cJ → e+ e- φ, φ → K+ K-
decay_card_for_signal = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3333 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c0
    1.000 e+ e- phi PHSP;
    Enddecay

    Decay chi_c1
    1.000 e+ e- phi PHSP;
    Enddecay

    Decay chi_c2
    1.000 e+ e- phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the combined signal (all three χ_cJ states)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3686_chicJ_ee_phi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_for_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "ChicJToeePhi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})       # ECMS = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# All three χ_cJ share the same final state (γ e+ e- K+ K-) and selection -> single Algorithm
event_selection = Selection.new
  .select_track {                     # Charged track selection: exactly four tracks
    cos_theta 0.93                    # |cosθ| < 0.93
    Vz        10.0                    # |Vz| < 10 cm
    Vr        1.0                     # Vr < 1 cm
    nChrp     "==2"                   # 2 positive tracks (e+, K+)
    nChrn     "==2"                   # 2 negative tracks (e-, K-)
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # Photon selection
    tdc_emc_start     0               # EMC TDC window [0, 14]
    tdc_emc_end       14
    angle_to_track    20.0            # angle to any charged track > 20 degrees
    energyThreshold_b 0.025           # barrel E > 25 MeV
    energyThreshold_e 0.050           # endcap E > 50 MeV
    nGam              ">=1"           # at least one photon
  }
  .pid(method: :probability) {        # PID by the probability method
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # high-p track -> lepton; EMC eraw>0.6 -> electron
    identify :kaon, against: [:pion, :proton]   # identify K+ and K-
    nlp "==1"                         # exactly one lepton+
    nlm "==1"                         # exactly one lepton-
    nkp "==1"                         # exactly one K+
    nkm "==1"                         # exactly one K-
  }
  # 4C kinematic fit of e+ e- K+ K- γ with four-momentum constraint
  .kinematic_fit([:gamma, :lp, :lm, :kp, :km]) {
    nominal                           # mark as the nominal fit
    constrain_four_momentum           # 4C energy-momentum constraint
    chi2_cut 200                      # loose χ² < 200 in BOSS (tight cut applied in ROOT)
  }

# Inexpressible BOSS-side procedures / assumptions
my_algorithm
  .note(:gamma_conversion_veto, "gamma-conversion veto applied before the 4C fit: e+e- pairs from photon conversion with vertex Rxy < 2 cm are rejected")
  .note(:pid_correction_method, "electron candidates required to satisfy 0.8 < E/p < 1.1 (EMC energy over MDC momentum); not expressible in the PID block")
  .note(:efficiency_curve, "photon energy required between 25 MeV and 1.4 GeV; the 1.4 GeV upper bound is applied on top of the 25 MeV threshold")
  .note(:kinematic_fit_tight_chi2, "nominal 4C fit uses a loose chi2 < 200 in BOSS; the paper's tight chi2 < 40 is applied later in ROOT")
  .with_decay_card(decay_card_for_signal).apply(event_selection)

# Execute on real data, inclusive MC, and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])