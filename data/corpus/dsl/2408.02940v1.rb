# DSL for BESIII paper 2408.02940v1: First observation of eta_c(2S) -> K+ K- eta
# psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K+ K- eta, eta -> gamma gamma
# Data: 2.712B psi(3686) events + continuum at sqrt(s)=3.650 GeV

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
# Continuum data at 3.650 GeV
continuum_3650 = DatasetManager.real_data.find("708_3650")
continuum_3650_incMC = DatasetManager.inclusive_mc.find("708_3650")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c(2S)              PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000 K+ K- eta                    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "etac2S_to_KK_eta"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("EtaC2S_KK_Eta")
alg.set_header(["EtaC2S_KK_EtaAlg/EtaC2S_KK_Eta.h"])
  .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz 10.0
                  Vr 1.0
                  nChrp ">=1"
                  nChrn ">=1"
                  nNet "==0"
                  nTot "==2"
                }
               .select_photon {
                  energyThreshold_b 0.025
                  energyThreshold_e 0.040  # Low threshold for M1 photon
                  angle_to_track 10.0
                  tdc_emc_start 0
                  tdc_emc_end 14
                  nGam ">=3"
                  nGam "<=6"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :kaon, against: [:pion, :proton]
                  nkp "==1"
                  nkm "==1"
                }
               # 4C kinematic fit (4-momentum conservation)
               .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma]) {
                  constrain_four_momentum
                  chi2_cut 20
                }
               # Competing 4-gamma hypothesis (to veto psi(3686) -> gamma gamma K+ K- eta)
               .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
                  constrain_four_momentum
                }
               # 5C kinematic fit with eta mass constraint
               .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 200
                }

alg.note(:eta_mass_window, "0.51 < M(gamma gamma) < 0.57 GeV/c^2 for eta signal region")
   .note(:radiative_photon, "photon with least energy selected as radiative photon from M1 transition")
   .note(:m4c_fit, "modified 4C: same as 5C but radiative photon energy allowed to float (suppress psi(3686)->K+K-eta with fake photon)")
   .note(:chi2_4c_vs_4gamma, "chi2_4c_3gamma < chi2_4c_4gamma to veto psi(3686)->gamma gamma K+K-eta")
   .note(:pi0_veto, "2D veto: |M(gamma_M1 gamma_low) - M(pi0)| < 0.025 GeV and |M(gamma_M1 gamma_high) - M(pi0)| < 0.04 GeV")
   .note(:etaprime_gamma_veto, "E(gamma2) outside [0.156, 0.196] GeV to suppress chi_c1 backgrounds")
   .note(:M3gamma_cut, "M(3gamma) > 0.6 GeV/c^2 to suppress K+K-pi0 backgrounds")
   .note(:KK_invmass_cut, "M(K+K-) < 3.0 GeV/c^2 to suppress J/psi bg; outside [1.007,1.033] GeV/c^2 to suppress phi bg")
   .note(:momentum_cut, "track momentum <= 2 GeV/c")
   .with_decay_card(decay_card)
   .apply(event_selection)

root_files = alg.execute_on([psip_data, psip_incMC, exMC_signal, continuum_3650, continuum_3650_incMC])