# ============================================================
# BESIII First Observation of h1(1380) in J/psi -> eta' KKbar pi
# Paper: arXiv:1804.05536v2
# Two decay modes:
#   Mode I:  J/psi -> eta' K+ K- pi0, eta' -> pi+ pi- eta, eta -> gamma gamma, pi0 -> gamma gamma
#            h1(1380) -> K*(892)+ K- + c.c.
#   Mode II: J/psi -> eta' Ks0 K+- pi-+, eta' -> pi+ pi- eta, eta -> gamma gamma
#            h1(1380) -> K*(892) Kbar + c.c.
# sqrt(s) = 3.097 GeV, 1.31e9 J/psi events
# ============================================================

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ========================================================================
# Mode I: J/psi -> eta' K+ K- pi0, eta' -> pi+ pi- eta, eta -> gamma gamma, pi0 -> gamma gamma
# Final state: pip pim kp km pi0 eta (= 4 charged + 4 photons)
# Focus: h1(1380) -> K*(892)+ K- + c.c.
# ========================================================================

decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.000 etap K+ K- pi0 PHSP;
    Enddecay
    Decay etap
    1.000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_etap_KKpi0_exMC"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("JpsiToEtapKKPi0")
alg_modeI.set_header(["JpsiToEtapKKPi0Alg/JpsiToEtapKKPi0.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
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
    nGam ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :pi0, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

alg_modeI
  .note(:etap_selection, "eta' selected by |M(pi+pi-eta)-m_etap|<0.03 GeV; pi0/eta from photon chi2 pairing: chi2_pi0eta < chi2_pi0pi0 AND chi2_pi0eta < chi2_etaeta")
  .note(:pi0_eta_mass_window, "pi0: |M(gamma gamma)-m_pi0|<0.02 GeV; eta: |M(gamma gamma)-m_eta|<0.03 GeV")
  .note(:kstar_selection, "K*(892)+- candidates selected via |M(K pi0)-m_K*|<0.15 GeV; K*(892)+ and K*(892)- constrained to equal yield")
  .note(:h1_extraction, "h1(1380) signal from simultaneous unbinned ML fit to K*(892) Kbar invariant mass: S-wave BW * q phase space factor; mass 1423.2+-2.1 MeV, width 90.3+-9.8 MeV")
  .note(:background_model, "Background from inclusive MC kernel estimation; phi eta eta and phi f0(1710) considered; etap sideband used")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID: highest C.L. hypothesis from TOF+dE/dx")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC TDC [0,14]*50ns; angle to track >10 deg")
  .note(:isospin_violation, "Isospin symmetry violation observed: B(h1->K*(892)+ K- + c.c.) vs B(h1->K*(892)0 Kbar0 + c.c.) differ due to K/K* mass differences near threshold")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])


# ========================================================================
# Mode II: J/psi -> eta' Ks0 K+- pi-+, eta' -> pi+ pi- eta, eta -> gamma gamma, Ks0 -> pi+ pi-
# Final state: pip pim kp km pim/pip Ks0 eta (= 6 charged + 2 photons)
# Focus: h1(1380) -> K*(892) Kbar + c.c. -> K*(892)+- K-+ and K*(892)0/barK*(892)0 Ks0/K+-
# ========================================================================

decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.000 etap K_S0 K+ pi- PHSP;
    Enddecay
    Decay etap
    1.000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_etap_KsKpi_exMC"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("JpsiToEtapKsKpi")
alg_modeII.set_header(["JpsiToEtapKsKpiAlg/JpsiToEtapKsKpi.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .with_decay_card(decay_card_modeII)

sel_modeII = Selection.new
sel_modeII.select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=3"
    nChrn ">=3"
    nNet "==0"
    nTot "==6"
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
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :K_S0, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

alg_modeII
  .note(:etap_selection, "eta' selected by |M(pi+pi-eta)-m_etap|<0.03 GeV; pion pair giving mass closest to etap nominal mass chosen")
  .note(:ks0_selection, "K_S0->pi+pi-: secondary vertex fit chi2<100; |M(pi+pi-)-m_KS0|<0.01 GeV; via min|M(pi+pi-)-m_KS0| if multiple candidates")
  .note(:eta_mass_window, "eta: |M(gamma gamma)-m_eta|<0.03 GeV")
  .note(:charged_pid, "3 pions and 1 kaon required from PID for the 4 non-KS0 tracks")
  .note(:kstar_selection, "K*(892) candidates: |M(Ks0 pi+-)-m_K*|<0.15 GeV, |M(K+- pi-+)-m_K*0|<0.15 GeV")
  .note(:h1_extraction, "h1(1380) signal from simultaneous unbinned ML fit to K*(892) Kbar invariant mass (combined Mode I + Mode II); mass 1423.2+-2.1 MeV, width 90.3+-9.8 MeV")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; KS0 daughters exempt from the standard selection")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC TDC [0,14]*50ns; angle to track >10 deg")
  .note(:isospin_violation, "Isospin symmetry violation observed: B(h1->K*(892)+ K- + c.c.) vs B(h1->K*(892)0 Kbar0 + c.c.) differ; measured both modes with Ks0 and KKpi0")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])