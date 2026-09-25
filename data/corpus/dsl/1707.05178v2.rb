# Search for eta_c -> pi+pi-pi0 and eta(1405) -> f0(980)pi0 in psi(3686) radiative decays
# BESIII, using 448.1 x 10^6 psi(3686) events

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay cards ###

# Analysis A: psi(3686) -> gamma eta_c, eta_c -> pi+ pi- pi0 (phase space)
decay_card_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c                              PHSP;
  Enddecay

  Decay eta_c
  1.000 pi+ pi- pi0                              PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                              PHSP;
  Enddecay

  End
DECAYCARD

# Analysis B: psi(3686) -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi-
decay_card_eta1405 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta(1405)                          PHSP;
  Enddecay

  Decay eta(1405)
  1.000 f_0(980) pi0                             PHSP;
  Enddecay

  Decay f_0(980)
  1.000 pi+ pi-                                  PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                              PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples ###
exMC_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_etac_3pi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_etac
  config.cross_section   = :default
end

exMC_eta1405 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_eta1405_f0pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_eta1405
  config.cross_section   = :default
end

### Common event selection (shared between the two analyses) ###
def build_selection
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93       # |cos(theta)| < 0.93
        Vz        10.0       # |Vz| < 10 cm (along beam)
        Vr        1.0        # |Vr| < 1 cm (transverse)
        nChrp     "==1"      # exactly one positive track
        nChrn     "==1"      # exactly one negative track
        nNet      "==0"      # net charge = 0
     }
     .select_photon {
        tdc_emc_start     0        # 0 <= TDC
        tdc_emc_end       14       # TDC <= 14 (50 ns/count)
        angle_to_track    10.0     # >=10 deg away from any charged track
        energyThreshold_b 0.025    # E > 25 MeV in EMC barrel  (|cos(theta)| < 0.80)
        energyThreshold_e 0.050    # E > 50 MeV in EMC end-caps (0.86 < |cos(theta)| < 0.92)
        nGam              ">=3"    # at least 3 photons
     }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]  # both tracks identified as pions using dE/dx (MDC) + TOF
        npip ">=1"
        npim ">=1"
     }
     # Main 4C fit under the gamma gamma gamma pi+ pi- hypothesis (nominal)
     .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 20                # chi2_4C(3gamma) < 20
     }
     # Competing 2-photon hypothesis (background suppression: stored chi2 for ROOT-level veto)
     .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
        constrain_four_momentum
     }
     # Competing 4-photon hypothesis (background suppression)
     .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {
        constrain_four_momentum
     }
  sel
end

### Algorithm for Analysis A: eta_c -> pi+ pi- pi0 ###
alg_etac = Algorithm.new("PsipGammaEtacTo3pi")
alg_etac.set_header(["PsipGammaEtacTo3piAlg/PsipGammaEtacTo3pi.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

sel_etac = build_selection

alg_etac
  .note(:pi0_selection,
        "pi0 candidate selected by |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2; " \
        "if multiple gamma gamma pairs pass, keep the one with M(gamma gamma) closest to nominal pi0 mass.")
  .note(:eta_veto,
        "veto events where the invariant mass of any of the other two gamma gamma photon pairs " \
        "lies within |M(gamma gamma) - m_eta| < 0.02 GeV/c^2, to reject events with an eta in the final state.")
  .note(:omega_veto,
        "reject events with |M(gamma pi0) - m_omega| < 0.05 GeV/c^2 to suppress omega -> gamma pi0 background.")
  .note(:jpsi_veto,
        "reject events with |M(gamma pi+ pi-) - m_J/psi| < 0.02 GeV/c^2 to suppress " \
        "psi(3686) -> pi0 J/psi (J/psi -> gamma pi+ pi- or J/psi -> pi+ pi- pi0 with a missing photon).")
  .note(:chi2_hypothesis_selection,
        "require chi2_4C(gamma gamma gamma pi+ pi-) < chi2_4C(gamma gamma pi+ pi-) AND " \
        "chi2_4C(gamma gamma gamma pi+ pi-) < chi2_4C(gamma gamma gamma gamma pi+ pi-) " \
        "to suppress backgrounds with two or four photons in the final state.")
  .note(:helix_correction,
        "helix parameter correction applied to charged pion tracks before the 4C kinematic fit; " \
        "correction parameters obtained from psi(3686) -> pi+ pi- pi0 control sample.")
  .with_decay_card(decay_card_etac)
  .apply(sel_etac)

alg_etac.execute_on([psip_data, psip_incMC, exMC_etac])

### Algorithm for Analysis B: eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi- ###
alg_eta1405 = Algorithm.new("PsipGammaEta1405ToF0Pi0")
alg_eta1405.set_header(["PsipGammaEta1405ToF0Pi0Alg/PsipGammaEta1405ToF0Pi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_eta1405 = build_selection

alg_eta1405
  .note(:pi0_selection,
        "pi0 candidate selected by |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2; " \
        "if multiple gamma gamma pairs pass, keep the one closest to nominal pi0 mass.")
  .note(:eta_veto,
        "|M(gamma gamma) - m_eta| > 0.02 GeV/c^2 for the other two photon pairs.")
  .note(:omega_veto,
        "|M(gamma pi0) - m_omega| > 0.05 GeV/c^2 to reject omega -> gamma pi0 background.")
  .note(:jpsi_veto,
        "|M(gamma pi+ pi-) - m_J/psi| < 0.02 GeV/c^2 vetoed to suppress psi(3686) -> pi0 J/psi backgrounds.")
  .note(:chi2_hypothesis_selection,
        "require chi2_4C(3-gamma pi+ pi-) < chi2_4C(2-gamma pi+ pi-) AND " \
        "chi2_4C(3-gamma pi+ pi-) < chi2_4C(4-gamma pi+ pi-).")
  .note(:helix_correction,
        "helix parameter correction applied to charged pion tracks before the 4C kinematic fit.")
  .with_decay_card(decay_card_eta1405)
  .apply(sel_eta1405)

alg_eta1405.execute_on([psip_data, psip_incMC, exMC_eta1405])
