# 2508.18761v1: chi_cJ -> Lambda anti-Lambda eta' (J=0,1,2) via psi(3686) -> gamma chi_cJ
# Two eta' decay modes:
#   Mode I: eta' -> gamma pi+ pi-
#   Mode II: eta' -> eta pi+ pi-, eta -> gamma gamma
# Lambda -> p pi-, anti-Lambda -> anti-p pi+
# psi(2S) data (BOSS 709_3686), 2.712 x 10^9 events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 HELAMP 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    Enddecay

    Decay chi_c0
    1.000 Lambda0 anti-Lambda0 eta' PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- DIY;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 HELAMP 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    Enddecay

    Decay chi_c0
    1.000 Lambda0 anti-Lambda0 eta' PHSP;
    Enddecay

    Alias another_eta eta

    Decay eta'
    1.000 another_eta pi+ pi- PHSP;
    Enddecay

    Decay another_eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chi_cJ_LLetap_ModeI_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chi_cJ_LLetap_ModeII_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

### Event selection — Mode I: eta' -> gamma pi+ pi- (final: p anti-p pi+ pi- pi+ pi- 2gamma) ###
alg_modeI = Algorithm.new("ChiCJLLEtap_ModeI")
alg_modeI.set_header(["ChiCJLLEtap_ModeIAlg/ChiCJLLEtap_ModeI.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=3"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
  nprm ">=1"
end
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
# Lambda reconstruction: p pi- secondary vertex fit
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# anti-Lambda reconstruction: anti-p pi+ secondary vertex fit
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# 4C kinematic fit: Lambda anti-Lambda gamma pi+ pi-
# Lambda and anti-Lambda mass windows enforced via invariant_mass constraints
.kinematic_fit([:Lambda, :Lambda_bar, :gamma, :pip, :pim]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:prp, :pim).between(1.110, 1.122)
  invariant_mass_of(:prm, :pip).between(1.110, 1.122)
  chi2_cut 18
end
.note(:lambda_mass_window,
  "Lambda signal region: |M(p pi-) - m_Lambda| < 5 MeV/c^2, decay length > 0. " \
  "anti-Lambda: |M(anti-p pi+) - m_anti-Lambda| < 5 MeV/c^2.")
.note(:background_veto,
  "Sigma^0 veto: M(Lambda/anti-Lambda gamma) > 1.2 GeV/c^2. " \
  "J/psi veto: |RM(pi+pi-) - m_J/psi| > 10 MeV/c^2, |RM(gamma1 gamma2) - m_J/psi| > 26 MeV/c^2. " \
  "pi0 veto: |M(gamma1 gamma2) - m_pi0| > 14 MeV/c^2. " \
  "eta' signal region: [m_etap - 12, m_etap + 12] MeV/c^2.")
.note(:competing_hypothesis,
  "Competing-fit vetoes: M(gamma pi+pi-) selection; photon choice minimizes |M(gamma pi+pi-) - m_etap|. " \
  "Non-eta' background from eta' sidebands [m_etap-60,m_etap-36] & [m_etap+36,m_etap+60] MeV/c^2.")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

### Event selection — Mode II: eta' -> eta pi+ pi-, eta -> gamma gamma ###
# (final: p anti-p pi+ pi- pi+ pi- 3gamma)
alg_modeII = Algorithm.new("ChiCJLLEtap_ModeII")
alg_modeII.set_header(["ChiCJLLEtap_ModeIIAlg/ChiCJLLEtap_ModeII.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=3"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=3"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
  nprm ">=1"
end
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Reconstruct eta from gamma pairs (kalman fit)
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end
# 5C kinematic fit: 4C + eta mass constraint
.kinematic_fit([:Lambda, :Lambda_bar, :eta, :pip, :pim]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:prp, :pim).between(1.110, 1.122)
  invariant_mass_of(:prm, :pip).between(1.110, 1.122)
  chi2_cut 53
end
.note(:lambda_mass_window,
  "Lambda signal region: |M(p pi-) - m_Lambda| < 5 MeV/c^2, decay length > 0.")
.note(:background_veto_modeII,
  "Sigma^0 veto: M(Lambda/anti-Lambda gamma) > 1.2 GeV/c^2. " \
  "eta' signal region: [m_etap - 10, m_etap + 10] MeV/c^2. " \
  "eta' sidebands: [m_etap-50,m_etap-30] & [m_etap+30,m_etap+50] MeV/c^2.")
.note(:signal_extraction,
  "Simultaneous unbinned maximum-likelihood fit to M(Lambda Lambda gamma pi+pi-) (Mode I) " \
  "and M(Lambda Lambda eta pi+pi-) (Mode II) in eta' signal and sideband regions. " \
  "chi_cJ mass windows: [3.390,3.445] (chi_c0), [3.495,3.525] (chi_c1), [3.540,3.570] (chi_c2) GeV/c^2. " \
  "Non-eta' background shape from eta' sidebands shared between signal/sideband regions. " \
  "Three chi_cJ signal MC shapes + Chebyshev polynomials for peaking/smooth backgrounds.")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])