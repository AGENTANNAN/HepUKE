# Paper 2503.21413v2: First observation of Λc(2595)+ and Λc(2625)+ → Λc+π⁰π⁰
# Data: 368.48 pb⁻¹ at √s=4.918 and 4.951 GeV
# Tag side: Λ̄_c⁻ with 3 modes: p̄K⁺π⁻, p̄K_S⁰, Λ̄π⁻
# Signal side: partial reconstruction of π⁰π⁰ (4γ, E_γ < 150 MeV)
# TagAnalysis: ST tag + partial_rec signal side

# ============================================================
# Decay cards
# ============================================================
decay_card_2595 = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_2625 = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Datasets (two energy points)
# ============================================================
data_4918 = DatasetManager.real_data.find("713_4918")
data_4951 = DatasetManager.real_data.find("713_4951")

incMC_4918 = DatasetManager.inclusive_mc.find("713_4918")
incMC_4951 = DatasetManager.inclusive_mc.find("713_4951")

# ============================================================
# Λc(2595)+ → Λc+π⁰π⁰
# ============================================================
alg_2595 = TagAnalysis.new("Lambdac2595ToLambdacPi0Pi0", "00-00-01")
  .set_header(["EventModel/Event.h", "EvtRecEvent/EvtRecTrack.h"])

alg_2595.tag_side(:"Lambda_c+") do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPi
end

alg_2595.signal_side do |s|
  s.photons 4
  s.at_least true
end

alg_2595.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_2595.note(:partial_reconstruction,
  "Partial reconstruction: only 4 photons from π⁰π⁰ are reconstructed. " \
  "Λc+ on signal side is NOT fully reconstructed. E_γ < 150 MeV per photon. " \
  "Signal extracted via fit to M_recoil(Λ̄_c⁻) distribution.")
alg_2595.note(:two_energy_points,
  "Data at √s=4.918 GeV (208.1 pb⁻¹) and √s=4.951 GeV (160.4 pb⁻¹) combined. " \
  "Two energy points handled via separate execute_on calls.")
alg_2595.note(:tag_modes,
  "Λ̄_c⁻ tag modes: p̄K⁺π⁻ (48.5%), p̄K_S⁰ (49.9%), Λ̄π⁻ (38.5%) efficiencies. " \
  "K_S⁰→π⁺π⁻, Λ→pπ⁻. M_tag ∈ (2.27,2.30) GeV/c². " \
  "M_K_S⁰ ∈ (0.487,0.511), M_Λ ∈ (1.111,1.121) GeV/c².")
alg_2595.note(:charge_conjugate,
  "Both Λ̄_c⁻ tag + Λc*+ signal and Λc+ tag + Λc*- signal included.")
alg_2595.note(:peaking_background,
  "Peaking background from Λc*+→Λc+π⁺π⁻ subtracted via MC.")

alg_2595.with_decay_card(decay_card_2595)
  .execute_on([data_4918, data_4951, incMC_4918, incMC_4951])

# ============================================================
# Λc(2625)+ → Λc+π⁰π⁰
# ============================================================
alg_2625 = TagAnalysis.new("Lambdac2625ToLambdacPi0Pi0", "00-00-01")
  .set_header(["EventModel/Event.h", "EvtRecEvent/EvtRecTrack.h"])

alg_2625.tag_side(:"Lambda_c+") do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPi
end

alg_2625.signal_side do |s|
  s.photons 4
  s.at_least true
end

alg_2625.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_2625.note(:partial_reconstruction,
  "Partial reconstruction: only 4 photons from π⁰π⁰ are reconstructed. " \
  "Λc+ on signal side is NOT fully reconstructed. " \
  "Signal extracted via fit to M_recoil(Λ̄_c⁻) distribution.")
alg_2625.note(:two_energy_points,
  "Data at √s=4.918 GeV (208.1 pb⁻¹) and √s=4.951 GeV (160.4 pb⁻¹) combined.")
alg_2625.note(:tag_modes,
  "Λ̄_c⁻ tag modes: p̄K⁺π⁻ (46.6%), p̄K_S⁰ (50.0%), Λ̄π⁻ (38.3%) efficiencies.")
alg_2625.note(:charge_conjugate,
  "Both Λ̄_c⁻ tag + Λc*+ signal and Λc+ tag + Λc*- signal included.")
alg_2625.note(:peaking_background,
  "Peaking background from Λc*+→Λc+π⁺π⁻ subtracted via MC.")

alg_2625.with_decay_card(decay_card_2625)
  .execute_on([data_4918, data_4951, incMC_4918, incMC_4951])