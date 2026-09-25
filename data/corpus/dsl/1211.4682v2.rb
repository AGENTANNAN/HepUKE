# =============================================================================
# BESIII: J/psi(psi') -> Lambda Lbar pi0 / Lambda Lbar eta
# arXiv:1211.4682v2
# =============================================================================

### Datasets ###
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")

### Decay cards ###
decay_card_jpsi_LLbar_pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_jpsi_LLbar_eta = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 eta PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_psip_LLbar_pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_psip_LLbar_eta = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda0 anti-Lambda0 eta PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_jpsi_LLbar_pi0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_LambdaLambdabar_pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_jpsi_LLbar_pi0
  c.cross_section   = :default
end

exMC_jpsi_LLbar_eta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_LambdaLambdabar_eta"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_jpsi_LLbar_eta
  c.cross_section   = :default
end

exMC_psip_LLbar_pi0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_LambdaLambdabar_pi0"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_psip_LLbar_pi0
  c.cross_section   = :default
end

exMC_psip_LLbar_eta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_LambdaLambdabar_eta"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_psip_LLbar_eta
  c.cross_section   = :default
end

# =============================================================================
# Common selection function to build the Selection chain
# =============================================================================
def build_selection
  sel = Selection.new
  sel.select_track {
       cos_theta 0.93
       # No impact-parameter cuts (tracks from secondary Lambda vertex)
       Vz  100.0
       Vr  10.0
       nChrp ">=2"
       nChrn ">=2"
     }
     .select_photon {
       nGam ">=2"
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start 0
       tdc_emc_end   14
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp ">=1"
       nprm ">=1"
     }
     .assign({:chrgp => :pip, :chrgn => :pim})
     .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     # 4C fit to J/psi(psi') -> Lambda Lbar gamma gamma
     .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
       nominal
       constrain_four_momentum
       chi2_cut 200   # loose; tight cut applied downstream in ROOT
     }
  sel
end

# =============================================================================
# ALG 1: J/psi -> Lambda Lbar pi0
# =============================================================================
alg1 = Algorithm.new("JpsiLLbarPi0")
alg1.set_header(["JpsiLLbarPi0Alg/JpsiLLbarPi0.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
alg1.note(:proton_pt_cut,
          "require pT(p) and pT(pbar) > 0.2 GeV/c; softer tracks are removed since MC "\
          "fails to describe them")
    .note(:sigma_plus_veto,
          "for J/psi -> Lambda Lbar pi0, remove |M(p pi0) - 1189.0| < 10 MeV/c^2 to veto "\
          "J/psi -> Sigma+ pi- Lbar (Sigma+ -> p pi0)")
    .note(:sigma0_sigma0bar_veto,
          "require M(Lambda Lbar) < 2.8 GeV/c^2 to remove J/psi -> Sigma0 Sbar0 background")
alg1.with_decay_card(decay_card_jpsi_LLbar_pi0).apply(build_selection)
alg1.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_LLbar_pi0])

# =============================================================================
# ALG 2: J/psi -> Lambda Lbar eta
# =============================================================================
alg2 = Algorithm.new("JpsiLLbarEta")
alg2.set_header(["JpsiLLbarEtaAlg/JpsiLLbarEta.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
alg2.note(:proton_pt_cut,
          "require pT(p),pT(pbar) > 0.2 GeV/c")
    .note(:sigma0_sigma0bar_veto,
          "require M(Lambda Lbar) < 2.6 GeV/c^2 to reject J/psi -> Sigma0 Sbar0")
alg2.with_decay_card(decay_card_jpsi_LLbar_eta).apply(build_selection)
alg2.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_LLbar_eta])

# =============================================================================
# ALG 3: psi(2S) -> Lambda Lbar pi0
# =============================================================================
alg3 = Algorithm.new("PsipLLbarPi0")
alg3.set_header(["PsipLLbarPi0Alg/PsipLLbarPi0.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
alg3.note(:proton_pt_cut,
          "pT(p),pT(pbar) > 0.2 GeV/c")
    .note(:psip_pipi_recoil_veto,
          "require |M_recoil(pi+ pi-) - 3097| > 8 MeV/c^2 to remove psi' -> pi+pi- J/psi background")
    .note(:sigma0_sigma0bar_veto,
          "require M(Lambda Lbar) < 3.08 GeV/c^2 to reject psi' -> gamma gamma J/psi "\
          "(J/psi -> Lambda Lbar) and psi' -> Sigma0 Sbar0")
alg3.with_decay_card(decay_card_psip_LLbar_pi0).apply(build_selection)
alg3.execute_on([psip_data, psip_incMC, exMC_psip_LLbar_pi0])

# =============================================================================
# ALG 4: psi(2S) -> Lambda Lbar eta
# =============================================================================
alg4 = Algorithm.new("PsipLLbarEta")
alg4.set_header(["PsipLLbarEtaAlg/PsipLLbarEta.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
alg4.note(:proton_pt_cut,
          "pT(p),pT(pbar) > 0.2 GeV/c")
    .note(:psip_pipi_recoil_veto,
          "require |M_recoil(pi+ pi-) - 3097| > 8 MeV/c^2 to remove psi' -> pi+pi- J/psi background")
    .note(:sigma0_sigma0bar_veto,
          "require M(Lambda Lbar) < 3.08 GeV/c^2")
alg4.with_decay_card(decay_card_psip_LLbar_eta).apply(build_selection)
alg4.execute_on([psip_data, psip_incMC, exMC_psip_LLbar_eta])
