# ============================================================
# BESIII First Observation of psi(3686) -> eta' e+ e-
# Paper: arXiv:1803.09714
# sqrt(s) = 3.686 GeV, 448.1M psi(3686) events
# Two eta' decay modes:
#   Mode I:  eta' -> gamma pi+ pi-
#   Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma
# ============================================================

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ========================================================================
# Mode I: psi(3686) -> eta' e+ e-, eta' -> gamma pi+ pi-
# Final state: 1 photon + e+ e- pi+ pi- (1 photon + 4 charged tracks)
# ========================================================================

decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 etap e+ e- PHSP;
    Enddecay
    Decay etap
    1.000 gamma pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2S_to_etap_ee_modeI_exMC"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("PsipToEtapEE_ModeI")
alg_modeI.set_header(["PsipToEtapEE_ModeIAlg/PsipToEtapEE_ModeI.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .with_decay_card(decay_card_modeI)

sel_modeI = Selection.new
sel_modeI.select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
    nTot "==4"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"
    nlm ">=1"
  }
  .remove([:lp <= :chrgp])
  .remove([:lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :lp, :lm, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

alg_modeI
  .note(:electron_pid, "Electron+positron require highest PID C.L. for electron hypothesis; pi+pi- no PID required. E/p>0.8 on higher-momentum track in e+e- pair for pi+pi-J/psi suppression")
  .note(:rm_pipi_cut, "Recoil mass of pi+pi- RM(pi+pi-) < 2.9 GeV/c2 to suppress psi(3686)->pi+pi-Jpsi, J/psi->e+e-; removes ~99.8% of this background")
  .note(:gamma_conversion_veto, "e+e- pair vertex distance from IP delta_xy < 2 cm to veto gamma conversion events from psi(3686)->etap gamma; removes >97% of conversion bkg")
  .note(:two_photon_veto, "cos_theta(e+) < 0.8 AND cos_theta(e-) > -0.8 to suppress two-photon process e+e-->e+e-etap")
  .note(:vertex_fit, "Vertex constraint on 4 charged tracks (pip pim ep em) to ensure they originate from IP, applied before 4C kinematic fit")
  .note(:helix_correction, "Helix parameter correction applied to charged tracks before 4C kinematic fit for better data/MC consistency")
  .note(:signal_fit, "M(gamma pi+ pi-) unbinned ML fit [0.85,1.05] GeV: MC-convolved Gaussian signal + 2nd-order Chebychev bkg + fixed peaking bkg (gamma conv + 2-photon)")
  .note(:tff_model, "Monopole TFF F(q2)=1/(1-q2/Lambda2) with Lambda=3.773 GeV; alternative Lambdas 3.2 and 5.0 GeV for systematics")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; photon: E>25 MeV barrel, E>50 MeV endcap; EMC TDC <700ns; angle to track >10deg")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])


# ========================================================================
# Mode II: psi(3686) -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
# Final state: 2 photons + e+ e- pi+ pi- (2 photons + 4 charged tracks)
# ========================================================================

decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 etap e+ e- PHSP;
    Enddecay
    Decay etap
    1.000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2S_to_etap_ee_modeII_exMC"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("PsipToEtapEE_ModeII")
alg_modeII.set_header(["PsipToEtapEE_ModeIIAlg/PsipToEtapEE_ModeII.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .with_decay_card(decay_card_modeII)

sel_modeII = Selection.new
sel_modeII.select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
    nTot "==4"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"
    nlm ">=1"
  }
  .remove([:lp <= :chrgp])
  .remove([:lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :lp, :lm, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

alg_modeII
  .note(:electron_pid, "Electron+positron require highest PID C.L. for electron hypothesis; pi+pi- no PID required")
  .note(:rm_pipi_cut, "Recoil mass of pi+pi- RM(pi+pi-) < 3.2 GeV/c2 to suppress psi(3686)->eta J/psi, J/psi->e+e-")
  .note(:eta_mass_window, "M(gamma gamma) in [0.520, 0.575] GeV/c2 additional eta mass window requirement")
  .note(:gamma_conversion_veto, "e+e- pair vertex distance from IP delta_xy < 2 cm to veto gamma conversion from psi(3686)->etap gamma")
  .note(:two_photon_veto, "cos_theta(e+) < 0.8 AND cos_theta(e-) > -0.8 to suppress two-photon process")
  .note(:vertex_fit, "Vertex constraint on 4 charged tracks to ensure they originate from IP before 4C kinematic fit")
  .note(:helix_correction, "Helix parameter correction applied to charged tracks before 4C kinematic fit")
  .note(:signal_fit, "M(gamma gamma pi+ pi-) unbinned ML fit [0.85,1.05] GeV: MC-convolved Gaussian signal + exponential bkg + fixed peaking bkg")
  .note(:tff_model, "Monopole TFF F(q2)=1/(1-q2/Lambda2) with Lambda=3.773 GeV")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; photon: E>25 MeV barrel, E>50 MeV endcap; EMC TDC <700ns; angle to track >10deg")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])