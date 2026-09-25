### Dataset description ###
# Paper: Study of h_c(1^1P_1) meson via psi(2S) -> pi0 h_c decays at BESIII
# arXiv: 2204.09413v2
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(2S) -> pi0 h_c, h_c -> gamma eta_c
# Inclusive channel: only pi0 recoil mass
# Tagged channel: additionally require E1 photon from h_c -> gamma eta_c
decay_card_hc = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c VSP_PWAVE;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_hc_pi0_gamma_etac"
  config.related_dataset = psip_data
  config.events = 300000
  config.decay_card = decay_card_hc
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PsiP2Pi0Hc"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz 10.0
                  Vr 1.0
                  nTot ">=2"       # suppress background
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end 14
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam ">=2"       # at least 2 for pi0 -> gamma gamma
                }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200
                  npi0 ">=1"
                }
                .note(:e1_photon_tag, "Tagged channel: additional E1 photon 465 < E_gamma < 535 MeV required; photon must not form pi0 with any other photon; if multiple pi0 in signal region RM(pi0) [3.500,3.550] GeV/c2, keep the one with minimum 1C fit chi2")
                .note(:jpsi_veto, "Veto psi(2S)->J/psi pi+pi- and psi(2S)->J/psi pi0 pi0: require pi+pi- recoil mass outside M_J/psi +/- 4 MeV/c2; pi0pi0 recoil mass outside M_J/psi [-8,+38] MeV/c2")
                .note(:total_energy_cut, "0.6 < E_EMC_Total < 3.2 GeV to suppress Bhabha and continuum backgrounds")
                .note(:kinematic_fit, "pi0 recoil mass RM(pi0) computed from 4-momentum conservation; signal extracted via fit to RM(pi0) distribution in ROOT analysis (convolution of resolution function + Breit-Wigner + 5th-order Chebyshev for background)")

my_algorithm.with_decay_card(decay_card_hc).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_hc])