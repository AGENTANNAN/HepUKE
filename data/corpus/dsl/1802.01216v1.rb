# BOSS Ruby DSL for BESIII paper 1802.01216v1
# Cross section measurement of e+e- → K+K-J/ψ and KS0KS0J/ψ at √s = 4.189-4.600 GeV
# Three decay modes (J/ψ → l+l-):
#   Mode I:   e+e- → K+K-J/ψ  (4C kinematic fit)
#   Mode II:  e+e- → KS0KS0J/ψ (6C kinematic fit with KS0 mass constraints; KS0 → π+π-)
#   Mode III: e+e- → π+π-J/ψ (4C kinematic fit; for cross section ratio)
#
# Multi-energy scan: uses 14 energy points.
# ECMS is NOT set in set_constant — injected at runtime per energy point.

### Dataset description ###
# Representative dataset; this analysis runs at multiple energy points 4.189-4.600 GeV
data_4260 = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# Decay card I: e+e- → K+K-J/ψ, J/ψ → l+l-
decay_card_KKJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000  K+  K-  J/psi  PHSP;
    Enddecay

    Decay J/psi
    1.000  e+  e-  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card II: e+e- → KS0KS0J/ψ, KS0 → π+π-, J/ψ → l+l-
decay_card_KSKSJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000  K_S0  K_S0  J/psi  PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-  PHSP;
    Enddecay

    Decay J/psi
    1.000  e+  e-  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card III: e+e- → π+π-J/ψ, J/ψ → l+l- (for ratio measurement)
decay_card_pipiJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  J/psi  PHSP;
    Enddecay

    Decay J/psi
    1.000  e+  e-  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples
exMC_KKJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_KKJpsi"
  config.related_dataset = data_4260
  config.events = 100000
  config.decay_card = decay_card_KKJpsi
  config.cross_section = :default
end

exMC_KSKSJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_KSKSJpsi"
  config.related_dataset = data_4260
  config.events = 100000
  config.decay_card = decay_card_KSKSJpsi
  config.cross_section = :default
end

exMC_pipiJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_pipiJpsi"
  config.related_dataset = data_4260
  config.events = 100000
  config.decay_card = decay_card_pipiJpsi
  config.cross_section = :default
end

### Event selection - Mode I: e+e- → K+K-J/ψ (J/ψ → l+l-) ###
# Charged tracks: ≥ 2 positive + 2 negative, |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm
# Leptons and kaons identified. 4C kinematic fit, χ²/dof < 10.
# Radiative Bhabha veto: cos(opening_angle) < 0.98 for all oppositely charged pairs.

alg_KKJpsi = Algorithm.new("KKJpsi")
alg_KKJpsi.set_header(["KKJpsiAlg/KKJpsi.h"])

sel_KKJpsi = Selection.new
sel_KKJpsi.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet  "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    nlp ">=1"; nlm ">=1"
  end
  .kinematic_fit([:kp, :km, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_KKJpsi
  .note(:helix_correction, "helix parameter correction applied; kinematic fit efficiency difference from control sample e+e-→K+K-π+π-: 3.8%")
  .note(:pid_correction_method, "tracking and PID efficiency corrections: 1.0% per pion, 2.5% per kaon; lepton tracking 1.0% per lepton")
  .note(:background_veto, "radiative Bhabha veto: cos(opening_angle) < 0.98 for all oppositely charged track pairs")
  .note(:efficiency_curve, "J/ψ signal/sideband method: signal [3.084,3.116] GeV/c², sidebands [3.004,3.068]+[3.132,3.196] GeV/c²; KK substructure efficiency weighting 10%")
  .with_decay_card(decay_card_KKJpsi)
  .apply(sel_KKJpsi)

### Event selection - Mode II: e+e- → KS0KS0J/ψ (KS0 → π+π-, J/ψ → l+l-) ###
# Charged tracks: ≥ 3 positive + 3 negative
# KS0 reconstruction: secondary vertex fit, L/σ > 4, M(π+π-) in [471,524] MeV/c², χ²_vertex < 100
# Leptons identified. 6C kinematic fit (4C + 2 KS0 mass constraints), χ²/dof < 10.

alg_KSKSJpsi = Algorithm.new("KSKSJpsi")
alg_KSKSJpsi.set_header(["KSKSJpsiAlg/KSKSJpsi.h"])

sel_KSKSJpsi = Selection.new
sel_KSKSJpsi.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"
    nChrn ">=3"
    nNet  "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    nlp ">=1"; nlm ">=1"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) do
    nominal
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
    invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
    constrain_four_momentum
    chi2_cut 200
  end

alg_KSKSJpsi
  .note(:helix_correction, "kinematic fit efficiency difference from control sample e+e-→KS0KS0π+π-: 5.9%")
  .note(:pid_correction_method, "tracking and PID efficiency: 1.0% per pion; lepton tracking 1.0%; KS0 reconstruction 3.0% per KS0")
  .note(:background_veto, "radiative Bhabha veto: cos(opening_angle) < 0.98 for all oppositely charged track pairs; KS0 flight significance L/σ > 4; veto e+e-→ππψ(3686) background")
  .note(:efficiency_curve, "KS0 mass window [471,524] MeV/c²; KS0 vertex fit χ² < 100; J/ψ signal/sideband method")
  .with_decay_card(decay_card_KSKSJpsi)
  .apply(sel_KSKSJpsi)

### Event selection - Mode III: e+e- → π+π-J/ψ (J/ψ → l+l-, for cross section ratio) ###
# Charged tracks: ≥ 2 positive + 2 negative
# 4C kinematic fit, χ²/dof < 10.

alg_pipiJpsi = Algorithm.new("PiPiJpsiXS")
alg_pipiJpsi.set_header(["PiPiJpsiXSAlg/PiPiJpsiXS.h"])

sel_pipiJpsi = Selection.new
sel_pipiJpsi.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet  "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    npip ">=1"; npim ">=1"
    nlp ">=1"; nlm ">=1"
  end
  .kinematic_fit([:pip, :pim, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_pipiJpsi
  .note(:helix_correction, "kinematic fit efficiency difference from control sample e+e-→π+π-π+π-: 2.6%")
  .note(:pid_correction_method, "tracking and PID: 1.0% per pion; lepton tracking 1.0% per lepton")
  .note(:background_veto, "radiative Bhabha veto: cos(opening_angle) < 0.98 for all oppositely charged track pairs")
  .note(:efficiency_curve, "J/ψ signal/sideband method; Zc(3900) substructure efficiency difference 4.0%")
  .with_decay_card(decay_card_pipiJpsi)
  .apply(sel_pipiJpsi)

### Execute ###
# Multi-energy scan: selection applied identically at all energy points
alg_KKJpsi.execute_on([data_4260, incMC_4260, exMC_KKJpsi])
alg_KSKSJpsi.execute_on([data_4260, incMC_4260, exMC_KSKSJpsi])
alg_pipiJpsi.execute_on([data_4260, incMC_4260, exMC_pipiJpsi])