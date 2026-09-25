# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data (2712 M events, 3.686 GeV)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC sample

# Decay card for the signal process:
# psi(3686) -> pi0 h_c, h_c -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c          PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-            PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000 gamma gamma      PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_pipijpsi_ee"
  config.related_dataset = psip_data          # associated real dataset
  config.events          = 100_000            # 100k events
  config.decay_card      = decay_card_signal  # signal decay card
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name     = "Pi0hc"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})           # E_cm = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                              # charged track selection
        cos_theta  0.93        # |cos(theta)| < 0.93
        Vz         100.0       # |Vz| < 100 cm
        Vr         10.0        # Vr < 10 mm
        nChrp      "==2"       # exactly 4 charged tracks: 2 positive ...
        nChrn      "==2"       # ... and 2 negative
        nNet       "==0"       # net charge zero
    }
    .select_photon {                             # photon selection
        tdc_emc_start      0       # TDC window 0 ...
        tdc_emc_end        14      # ... to 14
        angle_to_track     10.0    # min angle to nearest charged track (degrees)
        energyThreshold_b  0.025   # 25 MeV in the EMC barrel
        energyThreshold_e  0.050   # 50 MeV in the EMC endcap
        nGam               ">=2"   # at least two photons
    }
    .pid(method: :probability) {                 # probability-method PID
        # tracks with p > 1.0 GeV/c are treated as leptons;
        # electron if EMC energy > 1.0 GeV, otherwise muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 1.0
        identify :pion, against: [:kaon]         # remaining (p < 1.0 GeV/c) tracks: pi+ / pi- vs K
        npip "==1"                               # one pi+
        npim "==1"                               # one pi-
        nlp  "==1"                               # one l+
        nlm  "==1"                               # one l-
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {    # 1C Kalman fit to reconstruct pi0 -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                              # chi2 < 25
        npi0 ">=1"                               # at least one pi0 candidate
    }
    # Nominal 5C kinematic fit to pi+ pi- l+ l- pi0
    # (4C four-momentum conservation + the pi0 mass constraint already imposed above)
    .kinematic_fit([:pip, :pim, :lp, :lm, :pi0]) {
        nominal                  # nominal fit: its corrected four-momenta are the ones kept
        constrain_four_momentum  # total four-momentum constrained to the CMS energy
        chi2_cut 200             # loose cut in BOSS; the tight chi2 < 15 is applied in the ROOT analysis
    }
    # Alternative (competing) hypothesis: pi+ pi- pi0.
    # No chi2_cut and no nominal -> only this chi2 is stored, the veto is applied in ROOT.
    .kinematic_fit([:pip, :pim, :pi0, :lp, :lm]) {
        constrain_four_momentum
    }

# Generate the complete algorithm for the process defined in the decay card
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the signal exclusive MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])