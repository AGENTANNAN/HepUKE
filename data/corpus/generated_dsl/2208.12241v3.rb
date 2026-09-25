# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")        # psi(2S) data at sqrt(s) = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # matching inclusive MC sample

# Decay card: psi(2S) -> e+ e- eta_c  (eta_c decays inclusively)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 e+ e- eta_c PHSP;
    Enddecay

    End
DECAYCARD

# Dedicated exclusive signal MC for psi(2S) -> e+ e- eta_c (inclusive eta_c)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_ee_etac_inclusive"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_etac_recoil')

### Event selection (BOSS) ###
alg_name = "EtaCRecoil"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # sqrt(s) = 3.686 GeV

event_selection = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nTot      ">=2"     # at least two charged tracks
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    energyThreshold_b 0.025   # 25 MeV in the barrel
    energyThreshold_e 0.050   # 50 MeV in the endcap
    tdc_emc_start     0       # EMC shower time lower edge
    tdc_emc_end       14      # EMC shower time upper edge (14 x 50 ns = 700 ns)
  }
  .pid(method: :probability) {
    prob_cut 0.001                               # L(e) > 0.001
    identify :electron, against: [:pion, :kaon]  # e+/e- identified against pi and K
    nep ">=1"                                    # at least one e+
    nem ">=1"                                    # at least one e-  (i.e. at least one e+e- pair)
  }
  .remove(:ep) { condition "three_momentum_of(:ep) > 0.8" }   # keep e+ with p < 0.8 GeV/c
  .remove(:em) { condition "three_momentum_of(:em) > 0.8" }   # keep e- with p < 0.8 GeV/c
  .partial_rec([1, 2])   # recoil-mass technique: reconstruct e+ (recID 1) and e- (recID 2);
                         # the eta_c is inferred from the e+e- recoil mass (fit performed in ROOT)

my_algorithm
  .note(:pid_correction_method, "electron/positron identification uses L(e) > 0.001 together with the "
        "ratio requirement L(e)/(L(e)+L(pi)+L(K)) > 0.8; only L(e) > 0.001 maps onto prob_cut, the "
        "probability-ratio cut is applied inside the BOSS selection code")
  .note(:background_veto, "psi(3686) -> pi+pi- J/psi suppressed by rejecting any oppositely charged track "
        "pair, reinterpreted as pions, whose recoil mass lies in [3.090, 3.104] GeV/c^2")
  .note(:background_veto_pi0, "pi0 veto: gamma e+ e- combinations with invariant mass in "
        "[0.115, 0.150] GeV/c^2 are removed")
  .note(:background_veto_eta, "eta veto: gamma e+ e- combinations with invariant mass in "
        "[0.505, 0.570] GeV/c^2 are removed")
  .note(:background_veto_conversion, "gamma-conversion background rejected by requiring the e+ e- vertex "
        "R_xy < 2 cm and the e+ e- opening angle below 40 degrees")
  .note(:inclusive_etac_decay, "the eta_c is generated with inclusive decays; the eta_c left undecayed in "
        "the decay card is handled by the inclusive decay generator")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])