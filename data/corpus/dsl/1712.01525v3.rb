### Dataset description ###
# J/psi decays: 1.31 x 10^9 events at sqrt(s) = 3.097 GeV
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# Signal: J/psi -> gamma eta', eta' -> gamma pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma  eta'          PHSP;
    Enddecay

    Decay eta'
    1.0000  gamma  pi+  pi-      PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_gamma_etap_gammapipi"
    config.related_dataset = data_jpsi
    config.events = 1000000
    config.decay_card = decay_card_signal
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtapGammaPiPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                    cos_theta 0.93
                    Vz        10.0
                    Vr        1.0
                    nChrp     "==1"
                    nChrn     "==1"
                    nNet      "==0"
                }
                .select_photon {
                    tdc_emc_start    0
                    tdc_emc_end      14
                    angle_to_track   10.0
                    energyThreshold_b 0.040   # 40 MeV in barrel (this analysis)
                    energyThreshold_e 0.050
                    nGam             ">=2"
                }
                # 4C kinematic fit: J/psi -> gamma gamma pi+ pi-, chi2_4C < 100
                .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
                    nominal
                    constrain_four_momentum
                    chi2_cut 100
                }

my_algorithm
  .note(:pi0_veto, "|M(gamma gamma) - m(pi0)| > 20 MeV/c^2 required to remove J/psi -> pi+ pi- pi0 / gamma pi+ pi- pi0 background.")
  .note(:etap_selection, "Combination of gamma pi+ pi- closest to nominal eta' mass is kept as the eta' candidate; |M(gamma pi+ pi-) - m(eta')| < 20 MeV/c^2 required.")
  .note(:helix_correction, "Helix parameter correction applied to charged tracks in MC to match data 4C chi^2 distribution (Ref [34]).")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([data_jpsi, incMC_jpsi, exMC_signal])
