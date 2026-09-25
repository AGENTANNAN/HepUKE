# Measurement of J/ψ → Ξ(1530)⁻ Ξ̄⁺ and evidence for radiative decay Ξ(1530)⁻ → γ Ξ⁻
# arXiv:1911.06669v3 — single tag (ST) and double tag (DT) methods
# (1310.6 ± 7.0) × 10⁶ J/ψ events

### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### ========================================
### Mode A: ST — Ξ̄⁺ tag via partial reconstruction, Ξ(1530)⁻ missing
### ========================================
decay_card_ST = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi*- anti-Xi- J2BB3;
  Enddecay

  Decay anti-Xi-
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

exMC_ST = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_jpsi_Xi1530_barXi_inclusive"
  config.related_dataset = jpsi_data
  config.events = 300_000
  config.decay_card = decay_card_ST
  config.cross_section = :default
end

alg_ST = Algorithm.new("JpsiXi1530BarXiST")
alg_ST.set_header(["JpsiXi1530BarXiSTAlg/JpsiXi1530BarXiST.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

event_selection_ST = Selection.new
event_selection_ST.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"      # at least 2 π⁺ (one from Λ̄, one from Ξ̄⁺)
  nChrn ">=1"      # at least 1 p̄ from Λ̄ decay
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  identify :kaon, against: [:pion, :proton]
  nprm ">=1"
  npip ">=2"
end
# Remove identified kaons and protons; assign remaining tracks as pions
.remove([:kp <= :chrgp, :km <= :chrgn, :prp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
# Reconstruct Λ̄ → p̄ π⁺ via secondary vertex fit
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Partial reconstruction: miss Ξ(1530)⁻ (recID 1), reconstruct Ξ̄⁺ (recID 2)
# Ξ̄⁺ reconstructed by combining Λ̄ with remaining π⁺ via best_combination_by_mass
.partial_miss([1]) do
  best_combination_by_mass :Xi_bar, 1.32171
  require_recoil_mass 1.44, 1.65
end

alg_ST.note(:mass_window_Lambda, "M(p̄π⁺) required within 5 MeV/c² of nominal Λ mass (1.115683 GeV/c²); best candidate closest to nominal mass retained")
   .note(:mass_window_Xi_bar, "M(Λ̄π⁺) required within 8 MeV/c² of nominal Ξ⁺ mass (1.32171 GeV/c²)")
   .note(:decay_length, "Λ̄ and Ξ̄⁺ decay lengths required > 0 cm")
   .note(:charge_conjugate, "charge-conjugate mode (Ξ⁻ tag, Ξ̄(1530)⁺ recoil) included in actual BESIII analysis; this algorithm covers the primary Ξ̄⁺ tag mode")
   .note(:helix_correction, "helix parameter correction applied to charged tracks")
   .note(:continuum_veto, "continuum data at √s = 3.08 GeV used to verify negligible QED background")

alg_ST.with_decay_card(decay_card_ST).apply(event_selection_ST)
alg_ST.execute_on([jpsi_data, jpsi_incMC, exMC_ST])

### ========================================
### Mode B: DT — Ξ̄⁺ tag + Ξ⁻ + γ, radiative Ξ(1530)⁻ → γ Ξ⁻
### ========================================
decay_card_DT = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi*- anti-Xi- J2BB3;
  Enddecay

  Decay Xi*-
  1.0000 gamma Xi- PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda0 pi- PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Xi-
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

exMC_DT = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_jpsi_Xi1530_barXi_radiative"
  config.related_dataset = jpsi_data
  config.events = 300_000
  config.decay_card = decay_card_DT
  config.cross_section = :default
end

alg_DT = Algorithm.new("JpsiXi1530BarXiDT")
alg_DT.set_header(["JpsiXi1530BarXiDTAlg/JpsiXi1530BarXiDT.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

event_selection_DT = Selection.new
event_selection_DT.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"      # p + 2π⁺ (Λ̄ daughter + Ξ̄⁺ daughter)
  nChrn ">=3"      # p̄ + 2π⁻ (Λ daughter + Ξ⁻ daughter)
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  identify :kaon, against: [:pion, :proton]
  nprp ">=1"
  nprm ">=1"
  npip ">=2"
  npim ">=2"
end
.remove([:kp <= :chrgp, :km <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
# Reconstruct Λ̄ → p̄ π⁺ via secondary vertex fit (Ξ̄⁺ tag side)
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Reconstruct Λ → p π⁻ via secondary vertex fit (Ξ⁻ signal side)
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# 4C kinematic fit: γ + Ξ⁻(Λπ⁻) + Ξ̄⁺(Λ̄π⁺) → J/ψ
# Remaining pip/pim after second vertex fits are the Ξ̄⁺/Ξ⁻ daughter pions
.kinematic_fit([:gamma, :Lambda_bar, :pip, :Lambda, :pim]) do
  constrain_four_momentum
  chi2_cut 5
  nominal
end

alg_DT.note(:mass_window_Lambda, "M(pπ⁻) and M(p̄π⁺) required within 5 MeV/c² of nominal Λ mass (1.115683 GeV/c²); best candidate closest to nominal mass retained")
   .note(:mass_window_Xi, "M(Λπ⁻) and M(Λ̄π⁺) required within 8 MeV/c² of nominal Ξ mass (1.32171 GeV/c²); candidate with minimum |M - m_Ξ| selected")
   .note(:decay_length, "Λ, Λ̄, Ξ⁻, and Ξ̄⁺ decay lengths required > 0 cm")
   .note(:best_chisq_selection, "for each event, combination with lowest χ²_4C is selected")
   .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit")
   .note(:photon_barrel_endcap_gap, "photons in angular gap 0.80 < |cosθ| < 0.86 excluded (barrel/endcap transition)")
   .note(:photon_timing, "EMC timing 0 ≤ t ≤ 700 ns applied to suppress electronic noise")
   .note(:chi2_optimization, "χ²_4C < 5 determined by FOM = S/√(S+B) maximisation using signal and inclusive MC samples")
   .note(:charge_conjugate, "charge-conjugate mode included in actual BESIII analysis")

alg_DT.with_decay_card(decay_card_DT).apply(event_selection_DT)
alg_DT.execute_on([jpsi_data, jpsi_incMC, exMC_DT])