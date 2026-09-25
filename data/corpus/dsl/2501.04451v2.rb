# Paper 2501.04451v2 — Observation of the W-annihilation process D_s^+ → ωρ^+
# and measurement of D_s^+ → φρ^+ in D_s^+ → π^+π^+π^-π^0π^0 decays
# Double-tag method at √s = 4.128–4.226 GeV (7.33 fb⁻¹)
# Tag side: D_s^- → K_S^0 K^-, K^+ K^- π^-, K^+ K^- π^- π^0
# Signal side: D_s^+ → π^+ π^+ π^- π^0 π^0 (→ 4γ)

### Dataset preparation ###
# D_s energy scan data (BOSS 705 round12 + BOSS 703 round04/06/10)
ds_4130 = DatasetManager.real_data.find("705_4130")
ds_4160 = DatasetManager.real_data.find("705_4160")
ds_4180 = DatasetManager.real_data.find("703_4180")
ds_4190 = DatasetManager.real_data.find("703_4190")
ds_4200 = DatasetManager.real_data.find("703_4200")
ds_4210 = DatasetManager.real_data.find("703_4210")
ds_4220 = DatasetManager.real_data.find("703_4220")

ds_all = [ds_4130, ds_4160, ds_4180, ds_4190, ds_4200, ds_4210, ds_4220]

incMC_705 = DatasetManager.inclusive_mc.find("705_4130")
incMC_703 = DatasetManager.inclusive_mc.find("703_4180")

decay_card = <<~DECAYCARD
  Decay D_s+
  1.0000 pi+ pi+ pi- pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal exclusive MC for each energy point (D_s^+ → π^+π^+π^-π^0π^0)
exMC_signal = DatasetManager.create_exclusive_mc_for(ds_all) do |config|
  config.sample_name = "sig_Ds_pi2pi2pi0"
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (TagAnalysis — ST + signal side) ###
alg = TagAnalysis.new("DsToPiPiPiPi0Pi0")

alg.set_header(["DsToPiPiPiPi0Pi0Alg/DsToPiPiPiPi0Pi0.h"])
    .set_constant("ECMS" => [:double, 4.178])

alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0
  t.charm -1
end

alg.signal_side do |s|
  s.photons 4
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.min_photon_energy 0.025
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

alg.note(:multi_energy, "Data collected at √s = 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219 GeV; total integrated luminosity 7.33 fb⁻¹")
  .note(:tag_selection, "ST D_s^- tagged via 3 hadronic modes (K_S^0K^-, K^+K^-π^-, K^+K^-π^-π^0). Selection criteria for final-state particles (K_S^0, K^±, π^±, π^0, transition photon) and D_s^- candidates same as Ref. [50] (Phys. Rev. D 107, 052010)")
  .note(:kinematic_fit_6C, "6C kinematic fit: 4-momentum conservation + tag D_s^- mass → known D_s mass + (tag D_s + transition γ) mass → known D_s* mass. Multiple candidates per event: minimum χ²(6C) selected")
  .note(:kinematic_fit_7C, "7C kinematic fit: 6C + signal D_s^+ mass → known D_s mass. Updated four-momenta from 7C fit used for amplitude analysis")
  .note(:transition_photon, "Transition photon from D_s* → γ D_s: E_γ(lab) < 0.2 GeV required. Recoil mass M_rec in [1.95, 2.00] GeV/c² (5σ resolution window)")
  .note(:eta_veto, "η veto: events with any π^+π^-π^0 combination having invariant mass in [0.49, 0.58] GeV/c² (≈5σ resolution) are rejected; suppresses D_s^+ → π^+π^0η, η → π^+π^-π^0 background")
  .note(:Ks_veto_pi0pi0, "K_S^0 → π^0π^0 veto: M(π^0π^0) ∉ [0.487, 0.511] GeV/c² (5σ resolution)")
  .note(:Ks_veto_pipi, "K_S^0 → π^+π^- veto: secondary vertex fit on π^+π^- pairs; candidates with L/σ_L > 2 are rejected")
  .note(:D0_background_veto, "Open-charm D^0 background veto: events simultaneously satisfying |M(K^-π^+π^0) − M_D0| < 30 MeV/c² and |M(K^+π^+π^-π^-) − M_D0bar| < 30 MeV/c² are rejected (π/K misidentification). Analogous DD backgrounds excluded with same method")
  .note(:etap_veto, "η′ → π^+π^-γ background veto: two kinematic fits under signal (π^+π^+π^-π^0π^0) and background (ρ^+η′, η′→π^+π^-γ) hypotheses; event rejected if χ²(bkg) < χ²(sig)")
  .note(:signal_region, "D_s^+ invariant mass signal region [1.93, 1.99] GeV/c²; 1888 events retained with purity (79.3±1.3)%. Signal shape: MC shape ⊗ Gaussian; background shape from inclusive MC")
  .note(:amplitude_analysis, "Unbinned maximum-likelihood amplitude analysis with isobar formulation in covariant tensor formalism. 14 intermediate amplitudes retained (significance > 5σ each). Nonresonant component significance < 5σ")
  .note(:bf_measurement, "BF measured via DT method: B = Y_sig / Σ_i (Y_tag^i × ε_tag,sig^i / ε_tag^i). ST yields from Ref. [50]; DT yield 1985±68 from fit to M(D_s^+)")
  .with_decay_card(decay_card)
  .apply

root_files = alg.execute_on(ds_all + [incMC_705, incMC_703] + exMC_signal)