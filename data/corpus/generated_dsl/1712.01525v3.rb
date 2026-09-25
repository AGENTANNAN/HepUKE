### Dataset description ###
# J/psi peak real data (BOSS 7.0.8) and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> gamma eta', eta' -> gamma pi+ pi- (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 1M-event exclusive MC sample for the signal decay
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_gammapipi"
  config.related_dataset = jpsi_data      # associate with the real J/psi dataset
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaPrime"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})              # sqrt(s) = 3.097 GeV (J/psi peak)
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:helix_correction, "helix parameter corrections applied to MC charged tracks before the 4C kinematic fit to match data track resolution")

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {              # Charged track selection
                  cos_theta 0.93            # |cos(theta)| < 0.93
                  Vz        10.0            # |Vz| < 10 cm
                  Vr        1.0             # Vr < 1 cm
                  nChrp     "==1"           # exactly one positive track
                  nChrn     "==1"           # exactly one negative track
                  nNet      "==0"           # net charge zero
                }
               .select_photon {             # Photon selection
                  tdc_emc_start     0       # EMC TDC start
                  tdc_emc_end       14      # EMC TDC end
                  angle_to_track    10.0    # min angle to nearest charged track (degrees)
                  energyThreshold_b 0.040   # 40 MeV barrel threshold
                  energyThreshold_e 0.050   # 50 MeV endcap threshold
                  nGam              ">=2"   # at least two photons
                }
               .assign({:chrgp => :pip, :chrgn => :pim})  # no PID: hypotheses by charge, particle types fixed by the kinematic fit
               .kinematic_fit([:gamma, :gamma, :pip, :pim]) {   # 4C fit to gamma gamma pi+ pi-
                  nominal                                     # mark as the nominal fit
                  constrain_four_momentum                     # 4C energy-momentum constraint
                  chi2_cut 100                                # chi^2 < 100
                  invariant_mass_of(:gamma, :gamma).out_of(0.115, 0.155)    # pi0 veto: |M(gg) - m(pi0)| > 20 MeV
                  invariant_mass_of(:gamma, :pip, :pim).within(0.938, 0.978) # gamma pi+ pi- combination: |M - m(eta')| < 20 MeV
                }

# Generate and attach the algorithm for the decay card, then run on all datasets
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])