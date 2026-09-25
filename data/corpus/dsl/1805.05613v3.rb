# ============================================================
# BESIII Search for Invisible Decays of omega and phi
# Paper: arXiv:1805.05613
# Via J/psi -> V eta, eta -> pi+ pi- pi0
# sqrt(s) = 3.097 GeV, 1.31e9 J/psi events
# ============================================================

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ========================================================================
# Mode A: J/psi -> V(-> invisible) eta, eta -> pi+ pi- pi0
# Invisible decay search: V meson not reconstructed; signal appears in
# recoil mass against eta: RM(eta) = sqrt((Ecm - E_eta)^2 - P_eta^2)
# ========================================================================

decay_card_invisible = <<~DECAYCARD
    Decay J/psi
    1.000 omega eta PHSP;
    Enddecay
    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_inv = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_Veta_invisible_exMC"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_invisible
  config.cross_section = :default
end

alg_inv = Algorithm.new("JpsiToOmegaEtaInvisible")
alg_inv.set_header(["JpsiToOmegaEtaInvisibleAlg/JpsiToOmegaEtaInvisible.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .with_decay_card(decay_card_invisible)

sel_inv = Selection.new
sel_inv.select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=1"
    nChrn ">=1"
    nNet "==0"
    nTot "==2"
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
    identify :pion, against: [:kaon]
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_inv
  .note(:invisible_search, "Search for V(=omega,phi)->invisible via J/psi->V eta, eta->pi+pi-pi0; signal is a peak in M_recoil = sqrt((Ecm-E_eta)^2-P_eta^2) at V mass")
  .note(:eta_selection, "eta candidate: M(pi+pi-pi0) in [0.52,0.57] GeV; select combination with M(pi+pi-pi0) closest to eta nominal mass")
  .note(:vertex_fit, "Vertex fit on 2 charged pions to ensure common vertex origin")
  .note(:extra_photon_veto, "Sum of energies of extra photons (not used in eta reco) E_extra_gamma < 0.2 GeV to suppress V->gamma pi0, V->KsKL backgrounds")
  .note(:recoil_angle_cut, "Polar angle of system recoiling against eta: |cos_theta_recoil| < 0.7 to eliminate J/psi->X eta with X outside acceptance")
  .note(:fit_result, "Extended ML fit: signal PDF=MC shape, non-peaking=exponential, peaking fixed from MC (0.1 omega, 2.0 phi events). Nsig=1.4+-3.6(omega), -0.6+-4.5(phi) — consistent with zero")
  .note(:upper_limit, "90% CL upper limits: B(omega->inv)/B(omega->3pi) < 8.1e-5, B(phi->inv)/B(phi->KK) < 3.4e-4")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm, LED<10cm, LEr<1cm; PID: pion prob > kaon prob + electron prob")
  .note(:photon_selection, "E>25 MeV barrel(|cos_theta|<0.8), E>50 MeV endcap(0.86<|cos_theta|<0.92); EMC timing < 700ns; angle to track >10 deg")
  .note(:trigger_eff, "Trigger efficiency essentially 100% for events with 2 charged + 2 photons; 0.1% systematic assigned")
  .with_decay_card(decay_card_invisible)
  .apply(sel_inv)

alg_inv.execute_on([jpsi_data, jpsi_incMC, exMC_inv])


# ========================================================================
# Mode B (Reference): J/psi -> omega eta, omega -> pi+ pi- pi0, eta -> pi+ pi- pi0
# Visible reference channel for omega invisible decay normalization
# ========================================================================

decay_card_omega_vis = <<~DECAYCARD
    Decay J/psi
    1.000 omega eta PHSP;
    Enddecay
    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay
    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_omega_vis = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_omega_eta_visible_exMC"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_omega_vis
  config.cross_section = :default
end

alg_omega_vis = Algorithm.new("JpsiToOmegaEtaVisible")
alg_omega_vis.set_header(["JpsiToOmegaEtaVisibleAlg/JpsiToOmegaEtaVisible.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .with_decay_card(decay_card_omega_vis)

sel_omega_vis = Selection.new
sel_omega_vis.select_track {
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
    identify :pion, against: [:kaon]
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_omega_vis
  .note(:reference_channel, "Visible reference: J/psi->omega eta, omega->3pi, eta->3pi; used to normalize invisible decay upper limit")
  .note(:omega_eta_pairing, "2(pi+pi-pi0) final state: omega-eta pairing by chi2_omega_eta = (M_3pi_omega-m_omega)^2/sig_omega^2 + (M_3pi_eta-m_eta)^2/sig_eta^2; all 8 combinations tested, min chi2 selected")
  .note(:total_energy, "Selected candidate total energy E_tot > 2.95 GeV")
  .note(:recoil_angle_cut, "|cos_theta_recoil| < 0.7 for system recoiling against eta")
  .note(:pid_strategy, "PID for pions from eta decay only; pions from omega decay exempt (negligible omega->l+l-pi0 background)")
  .note(:mass_windows, "M(omega) in [0.65, 0.98] GeV; M(eta) in [0.41, 0.65] GeV")
  .note(:fit_2d, "2D ML fit: omega/eta signal = double CB; non-omega = 2nd-order Chebyshev; non-eta = reversed ARGUS. BKGIII (omega eta->gamma 3pi) = 1085.8+-126.6 evts, subtracted")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID: pion prob > kaon prob + electron prob")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC timing < 700ns; angle to track >10 deg")
  .with_decay_card(decay_card_omega_vis)
  .apply(sel_omega_vis)

alg_omega_vis.execute_on([jpsi_data, jpsi_incMC, exMC_omega_vis])


# ========================================================================
# Mode C (Reference): J/psi -> phi eta, phi -> K+ K-, eta -> pi+ pi- pi0
# Visible reference channel for phi invisible decay normalization
# ========================================================================

decay_card_phi_vis = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_phi_vis = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_to_phi_eta_visible_exMC"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_phi_vis
  config.cross_section = :default
end

alg_phi_vis = Algorithm.new("JpsiToPhiEtaVisible")
alg_phi_vis.set_header(["JpsiToPhiEtaVisibleAlg/JpsiToPhiEtaVisible.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .with_decay_card(decay_card_phi_vis)

sel_phi_vis = Selection.new
sel_phi_vis.select_track {
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
    identify :pion, against: [:kaon]
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .assign({:chrgp => :kp, :chrgn => :km})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_phi_vis
  .note(:reference_channel, "Visible reference: J/psi->phi eta, phi->K+K-, eta->3pi; used to normalize invisible decay upper limit")
  .note(:phi_selection, "phi candidate from 2 oppositely charged tracks assumed kaons without PID; M(K+K-) in [0.987, 1.10] GeV")
  .note(:eta_selection, "eta from pi+pi-pi0 with pi+pi- requiring pion PID, same procedure as invisible mode; M(3pi) closest to eta mass")
  .note(:total_energy, "Selected candidate total energy E_tot > 2.95 GeV")
  .note(:recoil_angle_cut, "|cos_theta_recoil| < 0.7 for system recoiling against eta")
  .note(:fit_2d, "2D ML fit: phi signal = rel. BW * Gaussian; eta signal = double CB; non-phi = reversed ARGUS; non-eta = reversed ARGUS. BKGIII (phi eta->gamma 3pi) = 238.6+-26.0 evts, subtracted")
  .note(:track_selection, "Charged: |cos_theta|<0.93, Vr<10mm, Vz<100mm; PID: pion prob > kaon prob + electron prob")
  .note(:photon_selection, "E>25 MeV barrel, E>50 MeV endcap; EMC timing < 700ns; angle to track >10 deg")
  .with_decay_card(decay_card_phi_vis)
  .apply(sel_phi_vis)

alg_phi_vis.execute_on([jpsi_data, jpsi_incMC, exMC_phi_vis])