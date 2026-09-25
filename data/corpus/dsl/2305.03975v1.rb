# Paper: 2305.03975v1
# Title: Measurement of CP-even fraction of D0 -> K_S0 pi+ pi- pi0
# Energy: 3.773 GeV (psi(3770), single energy point)
# Quantum-correlated DDbar pairs at psi(3770)
# Double-tag technique with CP-eigenstate and quasi-CP-eigenstate tag modes
# D0 -> K_S0 pi+ pi- pi0 (signal), tag D0bar in various modes

### Dataset preparation ###
data_712_3773 = DatasetManager.load_real_data.find("712_3773")
incMC_712_3773 = DatasetManager.load_inclusive_mc.find("712_3773")

all_data = [data_712_3773]
all_incMC = [incMC_712_3773]

# Decay card: psi(3770) -> D0 D0bar
# D0 -> K_S0 pi+ pi- pi0; D0bar -> K+ pi- (tag mode example)
# K_S0 -> pi+ pi-; pi0 -> gamma gamma
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0                           PHSP;
    Enddecay

    Decay D0
    1.000  K_S0  pi+  pi-  pi0                  PHSP;
    Enddecay

    Decay anti-D0
    1.000  K+  pi-                               PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                              PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                           PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "D0_KSpipipi0_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
# TagAnalysis: double-tag with quantum-correlated D0 D0bar pairs
alg = TagAnalysis.new("D0KSpipipi0_CP")
alg.set_header(["D0KSpipipi0CPAlg/D0KSpipipi0CP.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# D0bar tag side: use CP-eigenstate and quasi-CP tag modes
# Tag modes summarized in Table I of the paper
# CP-even tags: D0 -> K+K-, pi+pi-
# CP-odd tags: D0 -> K_S0 pi0, K_S0 eta, K_S0 omega, K_S0 eta'
# Quasi-CP tags: D0 -> K_S0 pi+ pi-, K_S0 pi+ pi- pi0, K_L0 pi0
d0_tag_modes = [
  :D0toKPi           # D0bar -> K+ pi- (flavor tag, not CP tag)
]

alg.tag_side(:D0) do |t|
  t.modes(*d0_tag_modes)
  t.charm -1            # anti-D0 tag
end

# Signal side: D0 -> K_S0 pi+ pi- pi0
alg.signal_side do |s|
  s.photons 2            # pi0 -> gamma gamma
  s.charged(pip: 3, pim: 1)  # K_S0 -> pi+ pi- + pi+ pi- from D0 decay
  s.require_charge 0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.with_decay_card(decay_card).apply

# Tag mode limitations in DTagAlg v1
alg.note(:tag_mode_unavailable,
  "Most CP-tag modes unavailable in DTagAlg v1: K+K-, pi+pi- (CP-even), K_S0 pi0, K_S0 eta, K_S0 omega, K_S0 eta' (CP-odd), K_S0 pi+ pi-, K_L0 pi0 (quasi-CP). Only D0->Kpi flavor tag available. Dropped from tag_side. Applied in ROOT.")

# Quantum-correlated DDbar analysis at psi(3770)
# CP-even fraction F_+ determined from double-tag yields
# Self-tag mode: D0 -> K_S0 pi+ pi- pi0 vs itself
# Signal: K_S0 -> pi+ pi- (secondary vertex), pi0 -> gamma gamma (mass-constrained fit)
# ST selection: |M(D)-M_D0| < 3sigma, DeltaE selection
# DT: M_BC and DeltaE requirements
alg.note(:quantum_correlation,
  "Quantum-correlated DDbar at psi(3770). CP-even fraction F_+ = 0.235+/-0.010_stat+/-0.002_syst from double-tag yields with CP-eigenstate and quasi-CP tag modes. Applied in ROOT.")

# D0 reconstruction: K_S0 + pi+ pi- + pi0
# K_S0: pi+pi- vertex fit, decay length > 2*sigma, M(pi+pi-) mass window
# pi0: gamma gamma with M(gammagamma) mass window and mass-constrained fit
# M_BC and DeltaE used for ST and DT selection
alg.note(:d0_reconstruction,
  "D0 reconstruction: K_S0 (pi+pi- vertex fit, L/sigma_L>2), pi+pi-, pi0 (gamma gamma mass-constrained fit). M_BC and DeltaE for ST/DT selection. Applied in ROOT.")

# ST yields from fits to M_BC distribution
# DT yields from mode-specific counting
# D0 D0bar mixing correction (1/(1-(2F_+^T-1)y_D)) applied to ST yields
alg.note(:yield_extraction,
  "ST yields from M_BC fits. DT from counting. D0-D0bar mixing correction applied. F_+ extracted from ratios of DT/ST yields. Applied in ROOT.")

# 2.93 fb^-1 at psi(3770)
alg.note(:luminosity,
  "2.93 fb^-1 at psi(3770). N_DDbar from previous BESIII measurement [11].")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)