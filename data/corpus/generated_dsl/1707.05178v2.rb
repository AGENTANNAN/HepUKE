# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card: psi(3686) -> gamma eta_c, eta_c -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_etac = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c        PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- pi0        PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi-, pi0 -> gamma gamma
decay_card_eta1405 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta(1405)    PHSP;
    Enddecay

    Decay eta(1405)
    1.0000 f0(980) pi0        PHSP;
    Enddecay

    Decay f0(980)
    1.0000 pi+ pi-            PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for each signal mode (500k events)
exMC_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gamma_etac_pip_pim_pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_etac
  config.cross_section   = :default
end

exMC_eta1405 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gamma_eta1405_f0_pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_eta1405
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Both signal modes (gamma eta_c and gamma eta(1405)) share the identical final state
# gamma gamma gamma pi+ pi- and the same selection criteria, so a single Algorithm /
# Selection chain serves both channels (shared final-state rule).
alg_name = "GammaEtaC"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93    # |cos(theta)| < 0.93
                  Vz          10.0    # |Vz| < 10 cm
                  Vr          1.0     # Vr < 1 cm in the transverse plane
                  nChrp       "==1"   # exactly one positive track
                  nChrn       "==1"   # exactly one negative track
                  nNet        "==0"   # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # TDC window 0-14 (50 ns/count)
                  tdc_emc_end       14
                  angle_to_track    10.0   # at least 10 deg from any charged track
                  energyThreshold_b 0.025  # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
                  energyThreshold_e 0.050  # E > 50 MeV in the endcaps (0.86 < |cos(theta)| < 0.92)
                  nGam              ">=3"  # at least three photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001                                    # probability > 0.001 (dE/dx + TOF based)
                  identify :pion, against: [:kaon, :proton]         # identify pi+ and pi- simultaneously
                  npip "==1"                                        # exactly one pi+
                  npim "==1"                                        # exactly one pi-
                }
               # Nominal 4C kinematic fit to gamma gamma gamma pi+ pi-
               .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
                  nominal                    # nominal fit: corrected four-momenta are saved
                  constrain_four_momentum    # 4C energy-momentum constraint
                  chi2_cut 20
                }
               # Competing hypothesis 2 gamma pi+ pi- : stores its chi2 for the ROOT-level veto
               .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
                  constrain_four_momentum
                }
               # Competing hypothesis 4 gamma pi+ pi- : stores its chi2 for the ROOT-level veto
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {
                  constrain_four_momentum
                }

# BOSS-side procedure that the DSL cannot express: helix-parameter correction of pion
# tracks (derived from a psi(3686) -> pi+ pi- pi0 control sample) applied before the fit
my_algorithm.note(:helix_correction,
  "helix parameter correction applied to all pion tracks, obtained from a psi(3686) -> pi+ pi- pi0 control sample, before the 4C kinematic fit")

# Attach one decay card (both cards describe the same final-state particles)
my_algorithm.with_decay_card(decay_card_etac).apply(event_selection)

# Execute on real data, inclusive MC, and both exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_etac, exMC_eta1405])