# BOSS DSL for: Lambda_c+ Weak Decay Asymmetry Measurements
# Paper: arXiv:1905.04707v1 (BESIII, 567 pb^-1 at 4.6 GeV)
# e+e- -> Lambda_c+ anti-Lambda_c-; partial reconstruction of one Lambda_c+
# Four signal modes: pK_S0, Lambda pi+, Sigma+ pi0, Sigma0 pi+
# This DSL expresses the pK_S0 mode (K_S0 -> pi+ pi-) as primary.
# Full angular analysis (helicity angles, simultaneous MLE fit over 4 modes)
# is performed at ROOT level.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4600 = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for e+e- -> Lambda_c+ anti-Lambda_c- at 4.6 GeV
# Lambda_c+ -> p K_S0 (signal), anti-Lambda_c- decays generically (unreconstructed)
# K_S0 -> pi+ pi-
decay_card_pKS = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 p+ K_S0 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p- K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_pKS = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "LambdacPKs_PartialRec"
  config.related_dataset = data_4600
  config.events = 500_000
  config.decay_card = decay_card_pKS
  config.cross_section = :default
end

# Algorithm: Lambda_c+ -> p K_S0, K_S0 -> pi+ pi- (primary mode)
alg_pKS = Algorithm.new("LambdacPKsAsym")
alg_pKS.set_header(["LambdacPKsAsymAlg/LambdacPKsAsym.h"])
        .set_constant({ "ECMS" => [:double, 4.600] })

sel_pKS = Selection.new

# Charged track selection: at least 3 tracks (p, pi+, pi- from K_S0)
# Paper: standard track quality |cos_theta|<0.93, |z|<10 cm, R<1 cm
# Lambda_c+ charge = +1, so net charge of (p, pi+, pi-) = +1
sel_pKS.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"    # p+ and pi+
  nChrn ">=1"    # pi-
  nNet "==1"     # net +1 = Lambda_c+ charge
end

# No photon requirement for pK_S0 mode
sel_pKS.select_photon do
  nGam ">=0"
end

# PID: identify proton, then remaining tracks are pions
# Paper uses dE/dx + TOF C.L. for pi/K/p hypotheses
sel_pKS.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:pion, :kaon]
  nprp ">=1"
end

sel_pKS.remove([:prp <= :chrgp])

# Remaining charged tracks (after proton removal) are pion candidates
sel_pKS.assign({ chrgp: :pip, chrgn: :pim })

# K_S0 -> pi+ pi- secondary vertex fit
sel_pKS.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Partial reconstruction: reconstruct Lambda_c+ from p + K_S0;
# anti-Lambda_c- is the missed (recoil) particle.
# M_BC signal window: [2.278, 2.294] GeV/c^2 (paper Fig.1)
# best_combination_by_mass selects candidate closest to nominal Lambda_c+ mass
# Decay card recID mapping: 0=psi(4260), 1=Lambda_c+, 2=anti-Lambda_c-,
#   3=p+, 4=K_S0, 5=anti-p-, 6=K_S0(bar side), 7=pi+, 8=pi-
# Miss anti-Lambda_c- (recID 2); all other valid recIDs are reconstructed
sel_pKS.partial_miss([2, 5, 6, 7, 8]) do
  best_combination_by_mass :Lambda_c, 2.28646
  require_recoil_mass 2.278, 2.294
end

# Inexpressible BOSS-side procedures captured as notes
alg_pKS.note(:partial_reconstruction, "Partial reconstruction method: only one Lambda_c+ is reconstructed; the anti-Lambda_c- partner is inferred from recoil mass (M_BC). The DSL uses partial_miss with the anti-Lambda_c- side as missed. Recoil mass window [2.278, 2.294] GeV/c^2 corresponds to the M_BC signal region.")
   .note(:best_candidate_selection, "If multiple Lambda_c+ candidates in an event, the one with smallest |Delta E| = |E_Lambda_c - E_beam| is kept. The DSL best_combination_by_mass uses mass proximity instead, which is an approximation of the Delta E ranking.")
   .note(:deltaE_requirement, "Delta E requirement is applied at ROOT level along with M_BC window. The exact window is not specified in the paper text but follows Ref. [25] conventions.")
   .note(:other_signal_modes, "Three additional Lambda_c+ decay modes studied: (1) Lambda pi+ with Lambda -> p pi-, (2) Sigma+ pi0 with Sigma+ -> p pi0 and pi0 -> gamma gamma, (3) Sigma0 pi+ with Sigma0 -> gamma Lambda and Lambda -> p pi-. Each mode requires its own Algorithm object with distinct track/photon/PID and vertex-fit chains. Only the pK_S0 mode is expressed in this DSL.")
   .note(:lambda_reconstruction, "Lambda -> p pi- reconstructed via secondary_vertex_fit, similar to K_S0 vertex fit. Required for Lambda pi+ and Sigma0 pi+ modes. Not expressed in this pK_S0 DSL.")
   .note(:pi0_reconstruction, "pi0 -> gamma gamma reconstructed via kalman_kinematic_fit with M(gamma gamma) constrained to nominal pi0 mass. Required for Sigma+ pi0 mode. Also uses pi0 mass window. Not expressed in this pK_S0 DSL.")
   .note(:sigma_plus_reconstruction, "Sigma+ -> p pi0 reconstructed from proton and pi0 candidates. Sigma+ has c*tau ~ 2.4 cm; vertex fit not typically used. Required for Sigma+ pi0 mode.")
   .note(:sigma_zero_reconstruction, "Sigma0 -> gamma Lambda, with Lambda -> p pi-. Sigma0 decays electromagnetically (c*tau ~ 2.2e-5 fm); the gamma and Lambda are combined without vertex fit. Required for Sigma0 pi+ mode.")
   .note(:ks_pi0pi0_veto, "K_S0 -> pi0 pi0 veto: M(pi0 pi0) outside [400, 550] MeV/c^2 applied in Sigma+ pi0 mode to suppress p K_S0 background where K_S0 -> pi0 pi0. Relevant only for Sigma+ pi0 mode.")
   .note(:angular_analysis, "Multi-dimensional simultaneous unbinned maximum likelihood fit over all four decay modes. Uses helicity formalism with angles theta_0, phi_1, theta_1, phi_2, theta_2, phi_3, theta_3. Free parameters: alpha_BP+, Delta_1^BP, sin Delta_0 (common). alpha_0 fixed to -0.20. Background subtraction via Type-I (wrong reconstruction) from inclusive MC and Type-II (combinatorial) from M_BC sideband [2.250, 2.270] GeV/c^2. MC-integration for normalization. All performed at ROOT level.")
   .note(:transverse_polarization, "First study of Lambda_c+ transverse polarization in unpolarized e+e- collisions. sin Delta_0 = -0.28 +/- 0.13(stat) +/- 0.03(syst), 2.1 sigma significance. PT(cos theta_0) = sqrt(1-alpha_0^2) cos theta_0 sin theta_0 sin Delta_0.")
   .note(:decay_asymmetry_results, "Measured alpha+ for: pK_S0 = 0.18 +/- 0.43(stat) +/- 0.14(syst) [first measurement]; Lambda pi+ = -0.80 +/- 0.11 +/- 0.02; Sigma+ pi0 = -0.57 +/- 0.10 +/- 0.07; Sigma0 pi+ = -0.73 +/- 0.17 +/- 0.07 [first measurement]. Sigma+ pi0 negative sign confirmed, contradicting positive model predictions by >8 sigma. Lambda pi+ and Sigma+ pi0 precision improved by factor ~3.")
   .note(:mBC_fit, "Signal yields from unbinned M_BC fits: signal shape + Type-I background (from signal MC, convolved with Gaussian for data-MC resolution difference) + Type-II background (ARGUS function). Signal region [2.278, 2.294], sideband [2.250, 2.270] GeV/c^2. ROOT-level procedure.")
   .note(:dataset_info, "Uses 567 pb^-1 at 4.6 GeV (BOSS 703, round07). DSL uses 703_4600 dataset.")
   .note(:cp_conservation, "CP conservation assumed: Delta_0 = Delta_0_bar, alpha_BP+ = -alpha_BbarPbar-, Delta_1^BP = -Delta_1^BP_bar. Input parameters: alpha_0 = -0.20 (from Ref.[27]), alpha_Lambda from BESIII (Ref.[22]), alpha_Sigma+ from PDG.")

alg_pKS.with_decay_card(decay_card_pKS).apply(sel_pKS)
alg_pKS.execute_on([data_4600, incMC_4600, exMC_pKS])