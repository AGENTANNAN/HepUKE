# frozen_string_literal: true
# =====================================================================
#  BOSS / BESIII DSL  --  confirmation of the charged charmoniumlike
#  Zc(3885)+-  in  e+e- -> pi+- (D D*)  via a double-D tag
#  (isospin processes I and II, with charge conjugates)
# =====================================================================

### =====================  Dataset preparation  =====================
# Real data and inclusive MC at the two energy points
#   4.230 GeV : 1092 pb^-1
#   4.260 GeV :  826 pb^-1
data_4230  = DatasetManager.real_data.find("703_4230")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
data_4260  = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

data_points = [data_4230, data_4260]   # both energy points

### =====================  Decay cards  =============================
# ---- Process I : e+e- -> pi+ D0 D*- , D*- -> anti-D0 pi- (soft pi- undetected)
decay_card_zc_proc1 = <<~DECAYCARD
  Decay psi(4260)
  1.0000  pi+  Zc(3885)-   PHSP;
  Enddecay

  Decay Zc(3885)-
  1.0000  D0  D*-   PHSP;
  Enddecay

  Decay D*-
  1.0000  anti-D0  pi-   PHSP;
  Enddecay

  End
DECAYCARD

# non-resonant phase-space card for process I
decay_card_phsp_proc1 = <<~DECAYCARD
  Decay psi(4260)
  1.0000  pi+  D0  anti-D0  pi-   PHSP;
  Enddecay

  End
DECAYCARD

# ---- Process II : e+e- -> pi+ D- D*0 , D*0 -> D0 pi0 (soft pi0 undetected)
decay_card_zc_proc2 = <<~DECAYCARD
  Decay psi(4260)
  1.0000  pi+  Zc(3885)-   PHSP;
  Enddecay

  Decay Zc(3885)-
  1.0000  D-  D*0   PHSP;
  Enddecay

  Decay D*0
  1.0000  D0  pi0   PHSP;
  Enddecay

  End
DECAYCARD

# non-resonant phase-space card for process II
decay_card_phsp_proc2 = <<~DECAYCARD
  Decay psi(4260)
  1.0000  pi+  D-  D0  pi0   PHSP;
  Enddecay

  End
DECAYCARD

# ---- D1(2420) Dbar0 background card
decay_card_d1bkg = <<~DECAYCARD
  Decay psi(4260)
  1.0000  pi+  D1(2420)-  anti-D0   PHSP;
  Enddecay

  Decay D1(2420)-
  1.0000  D*0  pi-   PHSP;
  Enddecay

  Decay D*0
  1.0000  D0  pi0   PHSP;
  Enddecay

  End
DECAYCARD

### =====================  Exclusive MC (200k events each)  ==========
# Zc signal MC and non-resonant phase-space MC for each isospin process,
# plus the D1(2420) Dbar0 background; one ExclusiveMC per energy point.
exMC_zc_proc1   = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "zc3885_sig_proc1"
  config.events        = 200_000
  config.decay_card    = decay_card_zc_proc1
  config.cross_section = :default
end

exMC_phsp_proc1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "zc3885_phsp_proc1"
  config.events        = 200_000
  config.decay_card    = decay_card_phsp_proc1
  config.cross_section = :default
end

exMC_zc_proc2   = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "zc3885_sig_proc2"
  config.events        = 200_000
  config.decay_card    = decay_card_zc_proc2
  config.cross_section = :default
end

exMC_phsp_proc2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "zc3885_phsp_proc2"
  config.events        = 200_000
  config.decay_card    = decay_card_phsp_proc2
  config.cross_section = :default
end

exMC_d1bkg      = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "d1_2420_d0bar_bkg"
  config.events        = 200_000
  config.decay_card    = decay_card_d1bkg
  config.cross_section = :default
end

all_exMC = [exMC_zc_proc1, exMC_phsp_proc1,
            exMC_zc_proc2, exMC_phsp_proc2,
            exMC_d1bkg].flatten

### ============================================================
### Process I : double-D tag  (D0 / anti-D0) + bachelor pi+
### ============================================================
alg_proc1 = TagAnalysis.new("Zc3885ProcI")
alg_proc1.set_header(["Zc3885ProcIAlg/Zc3885ProcI.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })   # measured CMS per run used in the 4C fit

# --- Tag side 1 : D0 / anti-D0 (both charges scanned)
alg_proc1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
end

# --- Tag side 2 : the recoiling D0 / anti-D0 of process I (double-D tag)
alg_proc1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
end

# --- Signal side : bachelor pi+ (net charge +1) not used by the D tags;
#     the soft pi- from the D*- is undetected.
alg_proc1.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge(1)
  s.min_photon_angle  20.0      # > 20 deg from any charged track
  s.min_photon_energy 0.025     # E > 25 MeV floor
  s.missing :pim                # soft pi- of D*- -> anti-D0 pi- is missing
end

# --- Kinematic fit : 4C with the D masses constrained to their PDG values
alg_proc1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200                # loose cut; the paper's tighter cut is applied downstream
end

alg_proc1
  .note(:track_selection, "all charged tracks required |cos(theta)| < 0.93 with PCA
    within 10 cm along the beam axis and 1 cm in the transverse plane; applied inside
    the DTagAlg tag reconstruction")
  .note(:pid_correction_method, "K/pi separation from combined dE/dx and TOF likelihoods:
    a track is called K if Prob(K) > Prob(pi); K_S0 daughters are exempt from this K/pi cut")
  .note(:pi0_reconstruction, "pi0 candidates built from photon pairs with
    0.115 < M(gamma gamma) < 0.150 GeV/c^2 and a 1C (mass) constraint")
  .note(:ks0_reconstruction, "K_S0 built from secondary-vertex-fitted pi+ pi- with
    0.487 < M(pi+ pi-) < 0.511 GeV/c^2 and a secondary-vertex constraint")
  .note(:signal_candidate_selection, "when several double-D tag combinations survive,
    keep the one with the smallest DeltaM_hat; require
    -20 < DeltaM_hat < 15 MeV/c^2 and |DeltaM| < 40 MeV/c^2 for the D0 D0bar tag
    (evaluated in ROOT on the stored tag observables)")
  .note(:background_veto, "veto D*D* events by |M_recoil(D pi) - M(D*)| < 30 MeV/c^2
    together with M(pi+ D0) >= 2.03 GeV/c^2 (applied downstream in ROOT)")
  .with_decay_card(decay_card_zc_proc1)
  .apply

### ============================================================
### Process II : double-D tag  (D- / D0) + bachelor pi+
### ============================================================
alg_proc2 = TagAnalysis.new("Zc3885ProcII")
alg_proc2.set_header(["Zc3885ProcIIAlg/Zc3885ProcII.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })

# --- Tag side 1 : D- (charge-conjugate D+, charm pinned to -1)
alg_proc2.tag_side(:Dplus) do |t|
  t.charm -1
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0,
          :DptoKsPiPiPi, :DptoKKPi
end

# --- Tag side 2 : the fully reconstructed D0 on the signal side (double-D tag)
alg_proc2.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
end

# --- Signal side : bachelor pi+ (net charge +1); the soft pi0 of the D*0 is undetected
alg_proc2.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge(1)
  s.photons 2                   # photons for pi0 reconstruction of the D0 modes
  s.min_photon_angle  20.0
  s.min_photon_energy 0.025
  s.missing :pi0                # soft pi0 of D*0 -> D0 pi0 is missing
end

# --- Kinematic fit : 4C with the D masses constrained to their PDG values
alg_proc2.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200                # loose cut; the paper's tighter cut is applied downstream
end

alg_proc2
  .note(:track_selection, "all charged tracks required |cos(theta)| < 0.93 with PCA
    within 10 cm along the beam axis and 1 cm in the transverse plane; applied inside
    the DTagAlg tag reconstruction")
  .note(:pid_correction_method, "K/pi separation from combined dE/dx and TOF likelihoods:
    a track is called K if Prob(K) > Prob(pi); K_S0 daughters are exempt from this K/pi cut")
  .note(:pi0_reconstruction, "pi0 candidates built from photon pairs with
    0.115 < M(gamma gamma) < 0.150 GeV/c^2 and a 1C (mass) constraint")
  .note(:ks0_reconstruction, "K_S0 built from secondary-vertex-fitted pi+ pi- with
    0.487 < M(pi+ pi-) < 0.511 GeV/c^2 and a secondary-vertex constraint")
  .note(:signal_candidate_selection, "when several double-D tag combinations survive,
    keep the one with the smallest DeltaM_hat; require
    -17 < DeltaM_hat < 14 MeV/c^2 and |DeltaM| < 35 MeV/c^2 for the D- D0 tag
    (evaluated in ROOT on the stored tag observables)")
  .note(:background_veto, "veto D*D* events by |M_recoil(D pi) - M(D*)| < 30 MeV/c^2
    together with M(pi+ D-) > 2.08 GeV/c^2 (applied downstream in ROOT)")
  .with_decay_card(decay_card_zc_proc2)
  .apply

### =====================  Execution  ================================
datasets = [data_4230, incMC_4230, data_4260, incMC_4260] + all_exMC

root_files_proc1 = alg_proc1.execute_on(datasets)
root_files_proc2 = alg_proc2.execute_on(datasets)