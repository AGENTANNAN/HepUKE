### Dataset preparation ###
# J/psi peak at 3.097 GeV: BOSS 7.0.8 sample with the full ~1.0e10 J/psi events.
jpsi_data  = DatasetManager.real_data.find("708_3097")       # Real J/psi data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Matching inclusive MC

# Decay card for J/psi -> gamma eta, eta -> pi+ pi- e+ e-  (EvtGen syntax)
decay_card_ee = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta       PHSP;
    Enddecay

    Decay eta
    1.000 pi+ pi- e+ e-   PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for J/psi -> gamma eta, eta -> pi+ pi- mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta         PHSP;
    Enddecay

    Decay eta
    1.000 pi+ pi- mu+ mu-   PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each of the two signal modes
exMC_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_gamma_eta_ee"
    config.related_dataset = jpsi_data          # tie simulation conditions to the real J/psi sample
    config.events          = 100_000
    config.decay_card      = decay_card_ee
    config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_gamma_eta_mumu"
    config.related_dataset = jpsi_data
    config.events          = 100_000
    config.decay_card      = decay_card_mumu
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaLL"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})      # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection chain for both lepton modes (identical final-state topology
# gamma pi+ pi- l+ l- and identical 4C fit hypothesis)
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta  0.93       # |cos(theta)| < 0.93
                  Vz         100.0      # |Vz| < 100 cm
                  Vr         10.0       # |Vr| < 10 mm in the transverse plane
                  nChrp      "==2"      # exactly two positive tracks
                  nChrn      "==2"      # exactly two negative tracks  -> four charged tracks
                  nNet       "==0"      # net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start      0      # EMC timing window start
                  tdc_emc_end        14     # EMC timing window end
                  angle_to_track     10.0   # angle to nearest charged track > 10 degrees
                  energyThreshold_b  0.025  # barrel energy threshold 25 MeV
                  energyThreshold_e  0.050  # endcap energy threshold 50 MeV
                  nGam               ">=1"  # at least one photon
                }
               .pid(method: :probability) {  # Particle identification (probability method)
                  prob_cut 0.001             # PID probability > 0.001
                  # high-momentum tracks (p > 1.0 GeV/c) treated as leptons;
                  # a lepton is called an electron if its EMC energy > 1.0 GeV, otherwise a muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 1.0
                  identify :pion, against: [:kaon]   # pion/kaon separation (pi+ and pi-)
                  npip "==1"                         # one pi+
                  npim "==1"                         # one pi-
                  nlp  "==1"                         # one l+ (e+ or mu+)
                  nlm  "==1"                         # one l- (e- or mu-)
                }
               # 4C kinematic fit to gamma pi+ pi- l+ l- (nominal fit).
               # Loose chi2 < 200 here; the tighter combined chi2_4C + sum(chi2_PID) < 50
               # is applied later in ROOT (Rule T3).
               .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
                  nominal                 # nominal fit: its corrected four-momenta are stored
                  constrain_four_momentum # 4C energy-momentum constraint to the CMS
                  chi2_cut 200
                }
               # Note: the photon-conversion veto (Phi_ee vs R_xy), the e/mu separation by
               # smallest combined chi2, the signal regions M(pi+pi-e+e-) in [0.53,0.57] GeV
               # and M(pi+pi-mu+mu-) in [0.531,0.567] GeV, and the ALP M(l+l-) selection are
               # post-kinematic-fit ROOT-level steps and are intentionally not encoded here.

# One algorithm serves both decay modes (same final-state topology and selection).
# The decay card fixes the kinematic-variable header definitions for the algorithm.
my_algorithm.with_decay_card(decay_card_ee).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_ee, exMC_mumu])