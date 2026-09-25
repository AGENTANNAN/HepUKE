### Dataset description ###
# XYZ data samples at 16 c.m. energy points from 4.009 to 4.600 GeV.
data_4260 = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# Signal decay card: e+e- -> pi0 pi0 psi(3686), psi(3686) -> pi+ pi- J/psi, J/psi -> l+ l-
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000  pi0  pi0  psi(2S)     PHSP;
    Enddecay

    Decay psi(2S)
    1.0000  pi+  pi-  J/psi       VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  e+   e-               PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000  gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000  pi0  pi0  psi(2S)     PHSP;
    Enddecay

    Decay psi(2S)
    1.0000  pi+  pi-  J/psi       VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  mu+  mu-              PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000  gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Background: e+e- -> pi+ pi- psi(3686), psi(3686) -> pi0 pi0 J/psi
decay_card_bkg = <<~DECAYCARD
    Decay psi(4260)
    1.0000  pi+  pi-  psi(2S)     PHSP;
    Enddecay

    Decay psi(2S)
    1.0000  pi0  pi0  J/psi       VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  e+   e-               PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000  gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_pi0pi0_psip_ee"
    config.related_dataset = data_4260
    config.events = 100000
    config.decay_card = decay_card_signal_ee
    config.cross_section = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_pi0pi0_psip_mumu"
    config.related_dataset = data_4260
    config.events = 100000
    config.decay_card = decay_card_signal_mumu
    config.cross_section = :default
end

exMC_bkg = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_pipi_psip_bkg"
    config.related_dataset = data_4260
    config.events = 100000
    config.decay_card = decay_card_bkg
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pi0pi0psip"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.416]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                    cos_theta   0.93
                    Vz          10.0
                    Vr          1.0
                    nChrp       "==2"
                    nChrn       "==2"
                    nNet        "==0"
                }
                .select_photon {
                    tdc_emc_start     0
                    tdc_emc_end       14
                    angle_to_track    10.0
                    energyThreshold_b 0.025
                    energyThreshold_e 0.050
                    nGam              ">=4"
                }
                .pid(method: :probability) {
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                   treat_as_electron_if_energy_above: 0.6
                    identify :pion, against: [:kaon]
                    npip "==1"
                    npim "==1"
                    nlp  "==1"
                    nlm  "==1"
                }
                # 4C kinematic fit: e+e- -> gamma gamma gamma gamma pi+ pi- l+ l-
                .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                    constrain_four_momentum
                    chi2_cut 120
                }
                # 7C kinematic fit: add mass constraints on two pi0 and J/psi
                .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                    nominal
                    constrain_four_momentum
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:'J/psi')
                    chi2_cut 200
                }

my_algorithm
  .note(:muon_hit_requirement, "Hit number requirement in the muon counter is applied on the mu+ mu- pair as in Ref [3]; not directly expressible in DSL.")
  .note(:photon_pairing, "Pairing of the four selected photons into the two pi0 candidates is chosen by minimizing (M(g1 g2)-M(pi0))^2 + (M(g3 g4)-M(pi0))^2.")
  .note(:pi0_mass_window, "|M(gamma_i gamma_j) - M(pi0)| < 20 MeV/c^2 required for each pi0 candidate.")
  .note(:jpsi_mass_window, "3.05 < M(l+ l-) < 3.15 GeV/c^2 required for J/psi candidates.")
  .with_decay_card(decay_card_signal_ee)
  .apply(event_selection)

root_files = my_algorithm.execute_on([data_4260, incMC_4260, exMC_signal_ee, exMC_signal_mumu, exMC_bkg])
