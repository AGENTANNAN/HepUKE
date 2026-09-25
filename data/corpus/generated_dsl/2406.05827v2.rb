### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data (three periods, ~20.3 fb^-1 in total)
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC sample

# Decay card for the Bhabha signal e+e- -> e+e- (phase space).
# No intermediate meson is produced, so psi(4260) is used as the top mother (KKMC/BESIII convention).
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.000 e+ e- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive Bhabha MC: 500k events, Babayaga@NLO cross section
exMC_bhabha = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_bhabha"
    config.related_dataset = psipp_data
    config.events          = 500_000
    config.decay_card      = decay_card_bhabha
    config.cross_section   = '/path/to/babayaga_nlo_xsec.dat'  # Babayaga@NLO cross section for e+e- -> e+e-
end

### Event selection (BOSS) ###
alg_name = "Bhabha"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})       # sqrt(s) = 3.773 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:ecm_calibration, "run-by-run calibrated E_cm applied for the three psi(3770) data periods
              (~20.3 fb^-1 in total); the integrated luminosity is measured from large-angle Bhabha scattering")
            .note(:electron_identification, "high-momentum tracks (p > 0.5 GeV/c) treated as leptons and
              electrons identified via E/p > 0.8; the DSL threshold maps the EMC-energy branch of this criterion")
            .note(:background_veto, "radiative Bhabha events with hard photons are vetoed using the defined
              photon candidates (25 MeV barrel / 50 MeV endcap, >10 deg from tracks); no explicit photon
              multiplicity requirement is applied in BOSS")

event_selection = Selection.new
event_selection.select_track {
                    cos_theta 0.93      # |cosθ| < 0.93
                    Vz        10.0      # |Vz| < 10 cm
                    Vr        1.0       # Vr < 1 cm in the transverse plane
                    nTot      "==2"     # exactly two charged tracks
                    nNet      "==0"     # net charge zero
                    nChrp     ">=1"     # at least one positive track
                    nChrn     ">=1"     # at least one negative track
                }
               .select_photon {
                    angle_to_track    10.0   # photon shower separated by >10 deg from the nearest charged track
                    energyThreshold_b 0.025  # 25 MeV in the EMC barrel
                    energyThreshold_e 0.050  # 50 MeV in the EMC endcap
                    # no nGam requirement here: photon candidates are only defined (hard-photon veto applied later)
                }
               .pid(method: :probability) {
                    prob_cut 0.001
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,   # p > 0.5 -> lepton
                                                   treat_as_electron_if_energy_above: 0.8    # electron identification
                    nlp ">=1"   # at least one positive lepton (e+)
                    nlm ">=1"   # at least one negative lepton (e-)
                }
               .kinematic_fit([:lp, :lm]) {
                    nominal                  # nominal fit: corrected four-momenta are saved
                    constrain_four_momentum  # 4C energy-momentum conservation on the e+e- hypothesis
                    chi2_cut 200             # loose cut; the optimal tight cut is applied in the ROOT analysis
                }

# Generate the complete algorithm for the process in the decay card
my_algorithm.with_decay_card(decay_card_bhabha).apply(event_selection)
# Execute on real data, inclusive MC and the exclusive Bhabha MC
root_files = my_algorithm.execute_on([psipp_data, psipp_incMC, exMC_bhabha])