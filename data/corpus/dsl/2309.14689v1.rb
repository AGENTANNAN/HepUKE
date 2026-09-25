# BESIII Analysis: psi(3686) → gamma eta_c(2S), eta_c(2S) → KKpi
# Paper: 2309.14689v1
# Data: psi(2S) at 3.686 GeV (BOSS 709) + continuum at 3.65 GeV (BOSS 709)

# ============================================================
# Datasets
# ============================================================
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")
cont_data   = DatasetManager.real_data.find("709_3650")
cont_incMC  = DatasetManager.inclusive_mc.find("709_3650")

# ============================================================
# Mode I: psi(3686) → gamma K_S0 K± pi∓, K_S0 → pi+ pi-
# Final state: gamma, pi+, pi-, K±, pi∓ (4 charged + 1 photon)
# ============================================================
decay_card_modeI = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma eta_c(2S)  VSP_PWAVE;
  Enddecay
  Decay eta_c(2S)
  1.0000 K_S0 K+ pi-  PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi2S_gamma_KSKpi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

alg_I = Algorithm.new("Etac2S_KSKpi_Analysis")
alg_I
  .set_header(["Etac2S_KSKpi_Analysis/Etac2S_KSKpi.h"])

sel_I = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
  end
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Modified 3C kinematic fit: fit with photon energy floating (to suppress fake-photon bg)
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_I
  .note(:ks_mass_window, "|M(pi+ pi-) - m_KS| < 7 MeV/c^2")
  .note(:ks_decay_length, "KS decay length > 2 * vertex resolution")
  .note(:modified_kinematic_fit, "Modified 3C fit: photon energy floats to suppress psi(3686)→KKpi background with fake photon. chi2_m3C < 20")
  .note(:ks_ks_veto, "Two non-KS charged tracks must have |M - m_KS| > 10 MeV to suppress gamma KS KS bg")
  .note(:recoil_mass_veto, "Recoil mass of pi+pi- pairs < 3.05 GeV to suppress J/psi backgrounds")
  .note(:fsr_correction, "FSR correction factor f_FSR applied to background shape from control samples")
  .note(:continuum_subtraction, "sqrt(s)=3.65 GeV continuum data used to estimate e+e-→gamma_ISR/FSR KKpi background")
  .note(:mass_fit, "Simultaneous extended unbinned maximum likelihood fit to M(KKpi) spectra for yield extraction")
  .with_decay_card(decay_card_modeI)
  .apply(sel_I)

alg_I.execute_on([psip_data, psip_incMC, signal_mc_I])

# ============================================================
# Mode II: psi(3686) → gamma K+ K- pi0, pi0 → gamma gamma
# Final state: gamma, K+, K-, gamma, gamma (2 charged + 3 photons)
# ============================================================
decay_card_modeII = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma eta_c(2S)  VSP_PWAVE;
  Enddecay
  Decay eta_c(2S)
  1.0000 K+ K- pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi2S_gamma_KKpi0"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

alg_II = Algorithm.new("Etac2S_KKpi0_Analysis")
alg_II
  .set_header(["Etac2S_KKpi0_Analysis/Etac2S_KKpi0.h"])

sel_II = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=3"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 15
    npi0 ">=1"
  end
  # Modified 4C kinematic fit: photon energy floats
  .kinematic_fit([:gamma, :kp, :km, :pi0]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_II
  .note(:modified_kinematic_fit, "Modified 4C fit: photon energy floats to suppress KKpi+gamma_fake bg. chi2_m4C < 15")
  .note(:pi0_mass_window, "pi0 mass window from kalman fit applied")
  .note(:omega_veto, "omega veto: |M_3gamma - m_omega| > 0.04 GeV to suppress psi(3686)→omega K+K- bg")
  .note(:muon_veto, "M_mumu < 2.90 GeV/c^2 to suppress J/psi→mu+mu- backgrounds")
  .note(:fsr_correction, "FSR correction from chi_c0 control sample")
  .note(:continuum_subtraction, "sqrt(s)=3.65 GeV data for continuum bg estimation")
  .note(:mass_fit, "Simultaneous fit to M(KKpi) with eta_c(2S) signal + chi_c1/2 + backgrounds")
  .with_decay_card(decay_card_modeII)
  .apply(sel_II)

alg_II.execute_on([psip_data, psip_incMC, signal_mc_II])