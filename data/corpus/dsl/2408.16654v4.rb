# Ξ⁰ → Λγ decay measurement at J/ψ
# Single-tag (ST) Ξ⁰→Λπ⁰ and double-tag (DT) Ξ⁰→Λγ
# J/ψ, 10087M events, arXiv:2408.16654v4

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
  Decay J/psi
  1.000 anti-Xi0 Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.000 anti-Lambda0 gamma PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "xi0_lambda_gamma"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("Xi0ToLambdaGamma")
alg.set_header(["Xi0ToLambdaGammaAlg/Xi0ToLambdaGamma.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

# ── Selection ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
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
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({ :chrgp => :pip, :chrgn => :pim })
  # Reconstruct Λ → pπ⁻ and Λbar → p̄π⁺ via secondary vertex fits
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct π⁰ → γγ via 1C kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Reconstruct ST Ξ⁰: choose the π⁰Λ combination with |M(π⁰Λ)-m_Ξ⁰| minimum
  # (minimum-mass-difference selection for Ξ⁰ candidate, performed in BOSS)
  .kinematic_fit([:Lambda, :pi0, :Lambda_bar, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 40
  }
  # Competing hypothesis for Σ⁰ background veto
  # (stored for ROOT-level |M(γΛbar)-M(Σ⁰)| > 12 MeV/c² cut)
  .kinematic_fit([:Lambda, :pi0, :Lambda_bar, :gamma]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

alg.with_decay_card(decay_card).apply(event_selection)

# ST yield from M_recoil fit (recoil mass against ST Ξ⁰bar); DT yield from p_Λ spectrum
# L/σ_L > 2 vertex separation, M(γΛbar)-M(Σ⁰) veto — applied in ROOT
alg.note(:st_dt_yield, "ST yield extracted from M_recoil fit recoiling against Ξ⁰bar; DT yield from p_Λ momentum spectrum in ROOT")
alg.note(:vertex_cut, "L/σ_L > 2 distance between primary vertex and Λ-Λbar intersection applied in ROOT")
alg.note(:sigma0_veto, "|M(γΛbar)-M(Σ⁰)| > 12 MeV/c² background veto applied in ROOT")

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])