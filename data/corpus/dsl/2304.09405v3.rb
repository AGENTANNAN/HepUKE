# Paper: 2304.09405v3
# Title: Measurement of branching fractions of Lambda_c+ decays to Sigma+ K+K-, Sigma+ phi, and Sigma+ K+ pi- (pi0)
# Energy: 7 energy points: 4.600-4.699 GeV
# Single-tag analysis: Lambda_c+ fully reconstructed; M_BC and DeltaE selection
# Signal modes relative to reference mode Lambda_c+ -> Sigma+ pi+ pi-
# Lambda_c+ -> Sigma+ K+K-, Sigma+ phi (phi -> K+K-), Sigma+ K+ pi-, Sigma+ K+ pi- pi0
# Sigma+ -> p pi0; pi0 -> gamma gamma

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")
incMC_703_4600 = DatasetManager.load_inclusive_mc.find("703_4600")

all_data = [data_703_4600]
all_incMC = [incMC_703_4600]

# Decay card: psi(4260) -> Lambda_c+ Lambda_c- (KKMC convention above threshold)
# Lambda_c+ -> Sigma+ K+K-; Sigma+ -> p pi0; pi0 -> gamma gamma
# Lambda_c- side not reconstructed (single tag)
decay_card_ref = <<~DECAYCARD
    Decay psi(4260)
    1.000  anti-Lambda_c-  Lambda_c+              PHSP;
    Enddecay

    Decay Lambda_c+
    1.000  Sigma+  pi+  pi-                       PHSP;
    Enddecay

    Decay Sigma+
    1.000  p+  pi0                                PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  pi-  pi+                      PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

exMC_ref = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_Sigma_pipi_reference"
  config.events         = 500000
  config.decay_card     = decay_card_ref
  config.cross_section  = :default
end

# Signal mode (example: Lambda_c+ -> Sigma+ K+K-)
decay_card_kk = <<~DECAYCARD
    Decay psi(4260)
    1.000  anti-Lambda_c-  Lambda_c+              PHSP;
    Enddecay

    Decay Lambda_c+
    1.000  Sigma+  K+  K-                         PHSP;
    Enddecay

    Decay Sigma+
    1.000  p+  pi0                                PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  pi-  pi+                      PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

exMC_kk = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_Sigma_KK_signal"
  config.events         = 500000
  config.decay_card     = decay_card_kk
  config.cross_section  = :default
end

### Event selection ###
# TagAnalysis: single-tag Lambda_c+ with pre-stored EvtRecDTag collection
# Reference mode: Lambda_c+ -> Sigma+ pi+ pi-
# Multiple signal modes each with their own TagAnalysis

# --- Reference mode: Lambda_c+ -> Sigma+ pi+ pi- ---

# Available Lambda_c+ tag modes. The paper reconstructs Lambda_c+ from scratch
# with M_BC and DeltaE, which is the ordinary Algorithm+Selection approach —
# not TagAnalysis. However, this IS a single-tag Lambda_c+ analysis.
# The Lambda_c+ modes used in DTagAlg are the hadronic Cabibbo-favored modes.

lambdac_modes = [
  :LambdacPtoPKPi,   # Lambda_c+ -> p K- pi+
  :LambdacPtoKsP,    # Lambda_c+ -> p K_S0
  :LambdacPtoKsPPiP, # Lambda_c+ -> p K_S0 pi+ pi-
  :LambdacPtoLambdPPiPPiM # Lambda_c+ -> Lambda pi+ pi- pi+
]

alg = TagAnalysis.new("Lc_Sigma_KK")
alg.set_header(["LcSigmaKKAlg/Lc_Sigma_KK.h"])
   .set_constant({"ECMS" => [:double, 4.600]})

alg.tag_side(:Lambdac) do |t|
  t.modes(*lambdac_modes)
  t.charm 1
end

# Signal side: Sigma+ K+ K-; Sigma+ -> p pi0, pi0 -> gamma gamma
# Or other signal modes
alg.signal_side do |s|
  s.photons 2            # pi0 -> gamma gamma
  s.charged(kp: 1, km: 1, prp: 1)  # K+ K- p (Sigma+ -> p pi0)
  s.require_charge 1
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.with_decay_card(decay_card_kk).apply

# This is a fully-reconstructed Lambda_c+ analysis using M_BC and DeltaE.
# The paper does NOT use the pre-stored tag collection (it reconstructs
# Lambda_c+ from scratch). The TagAnalysis framework provides the tag
# infrastructure but the paper's approach requires Algorithm+Selection.
# See tag-mode-vocabulary-gap memory.

alg.note(:full_reconstruction,
  "Paper reconstructs Lambda_c+ fully (not from pre-stored tags). M_BC and DeltaE used for signal extraction. The TagAnalysis tag_side declares modes for DTagAlg infrastructure; the actual reconstruction uses Algorithm+Selection approach. ROOT stage handles M_BC fits and DeltaE windows.")

# Multiple signal decay modes: Sigma+ K+K-, Sigma+ phi, Sigma+ K+ pi-, Sigma+ K+ pi- pi0
# Each mode analyzed separately; only one shown as example
alg.note(:signal_modes,
  "Five signal modes relative to reference Lambda_c+ -> Sigma+ pi+ pi-: (1) Sigma+ K+K-, (2) Sigma+ phi (phi->K+K-), (3) Sigma+ K+K- non-phi, (4) Sigma+ K+ pi-, (5) Sigma+ K+ pi- pi0. Each needs its own TagAnalysis or Algorithm. This file covers mode (1).")

# Sigma+ reconstruction: p pi0 with M(p pi0) in [1.174, 1.200] GeV/c^2 (~3 sigma)
# Lambda veto for pi+pi- modes: M(p pi-) outside [1.11, 1.12] GeV/c^2
# K_S0 veto: M(pi+pi-) outside [0.48, 0.52] GeV/c^2
alg.note(:sigma_selection,
  "Sigma+ -> p pi0: M(p pi0) in [1.174,1.200] GeV/c^2. Lambda veto: M(p pi-) outside [1.11,1.12]. K_S0 veto: M(pi+pi-) outside [0.48,0.52]. Applied in ROOT.")

# pi0: M(gammagamma) in [0.115,0.150] GeV/c^2; 1C kinematic fit constraining gamma gamma to pi0 mass; chi2 < 200
alg.note(:pi0_selection,
  "pi0: M(gammagamma) in [0.115,0.150] GeV/c^2; 1C kinematic fit with pi0 mass constraint, chi2 < 200. Applied in ROOT.")

# DeltaE requirements mode-dependent (Table 2):
# Sigma+ K+K-: -0.017 < DeltaE < 0.008
# Sigma+ K+ pi-: -0.014 < DeltaE < 0.008
# Sigma+ K+ pi- pi0: -0.028 < DeltaE < 0.012
# Sigma+ pi+ pi- (ref): -0.040 < DeltaE < 0.032
alg.note(:deltaE_selection,
  "Mode-dependent DeltaE requirements: Sigma+KK [-0.017,0.008]; Sigma+Kpi [-0.014,0.008]; Sigma+Kpipi0 [-0.028,0.012]; Sigma+pipi(ref) [-0.040,0.032] GeV. Candidate with minimal |DeltaE| accepted. Applied in ROOT.")

# M_BC signal extraction: unbinned ML fit with ARGUS background
# Signal: MC shape convolved with Gaussian
# 2D fit (M_BC vs M(K+K-)) for Sigma+ K+K- to separate phi and non-phi components
alg.note(:mbc_fit,
  "M_BC signal extraction: unbinned ML fit. Signal: MC shape conv. Gaussian; BG: ARGUS function. 2D fit (M_BC vs M_K+K-) for Sigma+KK mode to separate phi and non-phi. Simultaneous fit across 7 energy points with common RBF. Applied in ROOT.")

# Ratio of BF: RBF_ij = (N_i * epsilon_j * B_inter^j) / (N_j * epsilon_i * B_inter^i)
# Reference mode yields: 1123(47), 200(21), ... at 7 energy points
alg.note(:rbf_calculation,
  "RBF calculated as (N_sig * eps_ref * B_inter_ref) / (N_ref * eps_sig * B_inter_sig). Reference yields at 7 energy points from paper Table 3. Applied in ROOT.")

# Upper limit for Sigma+ K+ pi- pi0: likelihood scan with systematic uncertainties incorporated
alg.note(:upper_limit,
  "Upper limit on B(Lambda_c+ -> Sigma+ K+ pi- pi0) at 90% CL from likelihood scan incorporating systematic uncertainties. Applied in ROOT.")

# 7 energy points: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. 4.5 fb^-1 total. Signal yields in paper Table 3.")

# Charge-conjugate modes implied
alg.note(:charge_conjugate,
  "Throughout the analysis charge-conjugate modes are implicitly assumed.")

all_datasets = all_data + all_incMC + exMC_ref + exMC_kk
root_files = alg.execute_on(all_datasets)