# arXiv: 1901.03992v1
# Observation of X(3872) -> pi0 chi_cJ (J=0,1,2)
# BESIII Collaboration
# e+e- -> gamma X(3872) at 4.15-4.30 GeV, 9.0 fb^-1 total

# ============================================================================
# Dataset preparation
# ============================================================================

ds_4260 = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# ============================================================================
# Normalization channel: e+e- -> gamma X(3872) -> gamma pi+ pi- J/psi
# ============================================================================

decay_card_norm_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X_3872              PHSP;
  Enddecay

  Decay X_3872
  1.000  pi+ pi- J/psi             PHSP;
  Enddecay

  Decay J/psi
  1.000  e+ e-                     PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_norm_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X_3872              PHSP;
  Enddecay

  Decay X_3872
  1.000  pi+ pi- J/psi             PHSP;
  Enddecay

  Decay J/psi
  1.000  mu+ mu-                   PHOTOS VLL;
  Enddecay

  End
DECAYCARD

exMC_norm_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_X3872_pipiJpsi_ee"
  config.related_dataset = ds_4260
  config.events          = 500_000
  config.decay_card      = decay_card_norm_ee
  config.cross_section   = :default
end

exMC_norm_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_X3872_pipiJpsi_mumu"
  config.related_dataset = ds_4260
  config.events          = 500_000
  config.decay_card      = decay_card_norm_mumu
  config.cross_section   = :default
end

# ============================================================================
# Search channel: e+e- -> gamma1 X(3872) -> gamma1 pi0 chi_cJ, chi_cJ -> gamma2 J/psi
# ============================================================================

decay_card_search_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X_3872              PHSP;
  Enddecay

  Decay X_3872
  1.000  pi0 chi_c1                PHSP;
  Enddecay

  Decay chi_c1
  1.000  gamma J/psi               PHSP;
  Enddecay

  Decay J/psi
  1.000  e+ e-                     PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000  gamma gamma               PHSP;
  Enddecay

  End
DECAYCARD

decay_card_search_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X_3872              PHSP;
  Enddecay

  Decay X_3872
  1.000  pi0 chi_c1                PHSP;
  Enddecay

  Decay chi_c1
  1.000  gamma J/psi               PHSP;
  Enddecay

  Decay J/psi
  1.000  mu+ mu-                   PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000  gamma gamma               PHSP;
  Enddecay

  End
DECAYCARD

exMC_search_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_X3872_pi0chicJ_ee"
  config.related_dataset = ds_4260
  config.events          = 500_000
  config.decay_card      = decay_card_search_ee
  config.cross_section   = :default
end

exMC_search_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_X3872_pi0chicJ_mumu"
  config.related_dataset = ds_4260
  config.events          = 500_000
  config.decay_card      = decay_card_search_mumu
  config.cross_section   = :default
end

# ============================================================================
# Algorithm A: Normalization channel (gamma pi+ pi- l+ l-)
# No PID for pions; lepton identification via E/p
# 4C kinematic fit
# ============================================================================

alg_norm = Algorithm.new("X3872PiPiJpsi")
alg_norm.set_header(["X3872PiPiJpsiAlg/X3872PiPiJpsi.h"])
         .set_constant("ECMS" => [:double, 4.258])

sel_norm = Selection.new
sel_norm
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .select_photon {
    nGam ">=1"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_norm
  .note(:lepton_pid, "Electrons: E/p > 0.85; Muons: E/p < 0.25. DSL uses identify_high_momentum_leptons with p>1.0 and Eraw>0.6 thresholds as simplified proxy")
  .note(:chi2_optimized_cut, "Optimized chi2_4C/ndf < 10 from S/sqrt(S+B1+B2) maximization; chi2_cut 200 is a loose pre-cut")
  .note(:jpsi_mass_window, "M(l+l-) within +/-20 MeV/c^2 of J/psi nominal mass applied post-kinematic-fit; sideband regions 40 MeV/c^2 wide on each side with 20 MeV/c^2 gaps")
  .note(:mass_resolution_improvement, "M(pi+pi-J/psi) = M(pi+pi-l+l-) - M(l+l-) + M_0(J/psi) to improve resolution from 7.4 to 4.7 MeV/c^2")
  .note(:bhabha_veto, "cos(theta_pipi) < 0.98 and cos(theta_gamma_track) < 0.98 to suppress radiative Bhabha (e+e- -> e+e- gamma) background")
  .note(:eta_veto, "M(gamma pi+ pi-) > 0.6 GeV/c^2 and |M(gamma pi+ pi-) - M_0(eta')| > 0.02 GeV/c^2 to remove eta/eta' J/psi backgrounds")
  .note(:multiple_datasets, "Analysis combines 9.0 fb^-1 across ECM=4.15-4.30 GeV (signal region) plus 0.7 fb^-1 at 4.00-4.15 GeV and 2.8 fb^-1 at 4.30-4.60 GeV (sidebands)")
  .note(:isr_simulation, "ISR simulated with KKMC; cross-section lineshape from Ref.[12] used for signal MC weighting; FSR via PHOTOS")
  .note(:x3872_fit, "X(3872) yield from binned-likelihood fit to M(pi+pi-J/psi) with first-order polynomial background + MC-derived signal shape; significance from likelihood ratio")
  .with_decay_card(decay_card_norm_ee)
  .apply(sel_norm)

# ============================================================================
# Algorithm B: Search channel (gamma1 gamma2 pi0 l+ l-, pi0 -> gamma gamma)
# 4 photons + 2 leptons; 5C kinematic fit (4C + pi0 mass constraint)
# ============================================================================

alg_search = Algorithm.new("X3872Pi0ChicJ")
alg_search.set_header(["X3872Pi0ChicJAlg/X3872Pi0ChicJ.h"])
           .set_constant("ECMS" => [:double, 4.258])

sel_search = Selection.new
sel_search
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet "==0"
  }
  .select_photon {
    nGam ">=4"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }

alg_search
  .note(:lepton_pid, "Electrons: E/p > 0.85; Muons: E/p < 0.25")
  .note(:chi2_optimized_cut, "Optimized chi2_5C/ndf < 5 (4C + pi0 mass constraint); chi2_cut 200 is a loose pre-cut")
  .note(:jpsi_mass_window, "M(l+l-) within +/-20 MeV/c^2 of J/psi nominal mass applied post-kinematic-fit")
  .note(:pi0pi0_veto, "M(gamma1 gamma2) > 20 MeV/c^2 away from pi0 mass; same for gamma1/2 paired with higher-E pi0 daughter photon")
  .note(:omega_veto, "M(gamma1,2 pi0) < 0.732 GeV/c^2 to reject omega(782) -> gamma pi0 backgrounds")
  .note(:gamma_isr_veto, "Recoil mass against gamma1 or gamma2 both > 3.7 GeV/c^2 to suppress gamma_ISR psi(3686)")
  .note(:gamma_assignment, "gamma2 chosen to minimize |M(gamma2 J/psi) - M_0(chi_cJ)|; chi_cJ selection: DeltaM_0 < 25 MeV/c^2, DeltaM_{1,2} < 20 MeV/c^2")
  .note(:chi_cJ_loose_window, "Loose chi_cJ region: M(gamma1,2 J/psi) in [3.35, 3.60] GeV/c^2 applied for initial X(3872) search in M(pi0 chi_cJ)")
  .note(:mass_resolution_improvement, "M(pi0 chi_cJ) = M(pi0 l+l-) - M(l+l-) + M_0(J/psi) to improve mass resolution")
  .note(:gamma_interchange, "Signal MC includes events with gamma1/gamma2 interchanged; cross-feed among pi0 chi_cJ channels accounted for in fits")
  .note(:multiple_datasets, "Same multi-dataset combination as normalization channel")
  .note(:x3872_fit, "X(3872) yield from binned-likelihood fits to M(pi0 chi_cJ) for J=0,1,2 separately; signal shapes from MC including gamma1/2 interchange")
  .with_decay_card(decay_card_search_ee)
  .apply(sel_search)

# ============================================================================
# Execute
# ============================================================================

root_norm   = alg_norm.execute_on([ds_4260, incMC_4260, exMC_norm_ee, exMC_norm_mumu])
root_search = alg_search.execute_on([ds_4260, incMC_4260, exMC_search_ee, exMC_search_mumu])