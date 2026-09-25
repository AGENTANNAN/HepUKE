# e+e- → Ξ⁻ Ξ̄⁺ cross section measurement and Ξ(1820) observation
# arXiv:1910.04921v2 — single baryon tag method
# 11.0 fb⁻¹ at √s = 4.009–4.6 GeV, 15 energy points

### Dataset description ###
data_points = [
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4600"),
]

incMC_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

### Decay card for e+e- → Ξ⁻ Ξ̄⁺ (KKMC + psi(4260) top mother) ###
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Xi- anti-Xi- PHSP;
  Enddecay

  Decay Xi-
  1.0000 pi- Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Xi-
  1.0000 pi+ anti-Lambda0 PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

### Exclusive MC for energy scan ###
exMC = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_Xi_Xi_bar"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("XiXiBarAnalysis")
alg.set_header(["XiXiBarAnalysisAlg/XiXiBarAnalysis.h"])

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  nChrp ">=1"      # at least 1 proton from Λ decay
  nChrn ">=2"      # at least 2 π⁻ (one from Λ, one from Ξ⁻)
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  identify :kaon, against: [:pion, :proton]
  nprp ">=1"
  npim ">=2"
end
# Remove identified kaons and protons from generic charged lists
.remove([:kp <= :chrgp, :km <= :chrgn, :prp <= :chrgp, :prm <= :chrgn])
# Assign remaining charged tracks as pions
.assign({chrgp: :pip, chrgn: :pim})
# Reconstruct Λ → p π⁻ via secondary vertex fit
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Partial reconstruction: miss Ξ̄⁺ (recID 2), reconstruct Ξ⁻ from π⁻ + Λ
.partial_miss([2]) do
  best_combination_by_mass :Xi, 1.32171
  require_recoil_mass 1.2, 1.5
end

# Attach notes for BOSS-side procedures not expressible in formal DSL constructs
alg.note(:secondary_vertex_chi2, "Λ secondary vertex fit requires χ² < 500 (3 d.o.f.); Ξ⁻ secondary vertex fit performed via best_combination_by_mass with mass-difference minimisation")
   .note(:mass_window_Lambda, "M(pπ⁻) required within 5 MeV/c² of nominal Λ mass (1.115683 GeV/c²)")
   .note(:mass_window_Xi, "M(π⁻Λ) required within 10 MeV/c² of nominal Ξ⁻ mass (1.32171 GeV/c²)")
   .note(:decay_length, "Λ and Ξ⁻ decay lengths required > 0 cm")
   .note(:charge_conjugate, "charge-conjugate mode (Ξ̄⁺ tag, Ξ⁻ recoil) included in actual BESIII analysis; this algorithm covers the primary Ξ⁻ tag mode")

alg.with_decay_card(decay_card).apply(event_selection)
root_files = alg.execute_on(data_points + incMC_points + exMC)