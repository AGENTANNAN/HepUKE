# 1907.13041v2: ψ(3686) → Ξ(1530)- anti-Ξ(1530)+ and Ξ(1530)- anti-Ξ+
# Single-baryon tagging technique with partial reconstruction

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3686_data  = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → Ξ(1530)- anti-Ξ(1530)+ with full decay chain
# Tag side (reconstructed): Ξ(1530)- → π- Ξ0 → π- π0 Λ → π- π0 p π- → π- γγ p π-
# Recoil side (missed): anti-Ξ(1530)+ → π+ anti-Ξ0 → π+ π0 anti-Λ → π+ γγ anti-p π+
decay_card_Xi1530 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Xi(1530)- anti-Xi(1530)+ PHSP;
  Enddecay
  Decay Xi(1530)-
  1.0000 pi- Xi0 VSS;
  Enddecay
  Decay anti-Xi(1530)+
  1.0000 pi+ anti-Xi0 VSS;
  Enddecay
  Decay Xi0
  1.0000 pi0 Lambda0 PHSP;
  Enddecay
  Decay anti-Xi0
  1.0000 pi0 anti-Lambda0 PHSP;
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

exMC_Xi1530 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Xi1530_tag"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card_Xi1530
  config.cross_section = :default
end

# ====== Algorithm: single-baryon tag of Ξ(1530)-, partial recoil for anti-baryon ======
alg = Algorithm.new("Xi1530Tag")
alg.set_header(["Xi1530TagAlg/Xi1530Tag.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })

# Selection: reconstruct tag-side (Ξ(1530)- → π-Ξ0), infer anti-baryon from recoil
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=1"
  nChrn ">=2"
}
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  tdc_emc_start 0
  tdc_emc_end 700
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  npim ">=2"
}
.remove([:prp <= :chrgp])
.remove([:prm <= :chrgn])
.assign({:chrgp => :pip, :chrgn => :pim})
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 20
  npi0 ">=1"
}
.secondary_vertex_fit([:prp, :pim]) {
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
}
.partial_miss([2]) {
  best_combination_by_mass :Xi0, 1.31486
  best_combination_by_mass :"Xi(1530)-", 1.535
  require_recoil_mass 1.20, 1.70
}

alg.with_decay_card(decay_card_Xi1530).apply(event_selection)
alg.note(:tag_modes, "Two tag modes used in full analysis: π-Ξ0 and π0Ξ-; π0Ξ- mode requires separate decay card. Results from both modes combined.")
    .note(:cascades, "Ξ0 candidates from π0Λ within 10 MeV/c^2 of nominal Ξ mass; Ξ- candidates from π-Λ within 10 MeV/c^2 with secondary vertex fit")
    .note(:Lambda_selection, "Λ → pπ- within 5 MeV/c^2 of nominal mass; secondary vertex fit χ^2 < 500; positive decay length required")
    .note(:pi0_window, "π0 1C fit χ^2 < 20; mass window and further cuts in ROOT")
    .note(:Xi1530_window, "M(πΞ) within 15 MeV/c^2 of nominal Ξ(1530)- mass; best candidate per event")
    .note(:Xi_decay_length, "Ξ- decay length required positive in ROOT; multiple candidate resolution by closest mass")
    .note(:recoil_fit, "Signal yields from extended ML fit to M(πΞ)^recoil spectrum; signal shape from MC convolved with Gaussian; WCB from MC shape; Other-Bkg: 3rd-order Chebychev")
    .note(:alpha_measurement, "Angular distribution parameter α from cosθ_B fit in [-0.8, 0.8], 8 bins; efficiency-corrected yields per bin from M^recoil fits")
    .note(:double_counting, "~10% double-counting of Ξ(1530)-anti-Ξ(1530)+ final state from single-baryon method; accounted for in combination")
    .note(:off_peak, "Continuum background estimated from 3.65 GeV off-peak data (44 pb-1); found negligible")
    .note(:cc_mode, "Charge-conjugate mode (anti-Ξ(1530)+ tag) uses same algorithm with charge-conjugated selection")

alg.execute_on([psi3686_data, psi3686_incMC, exMC_Xi1530])