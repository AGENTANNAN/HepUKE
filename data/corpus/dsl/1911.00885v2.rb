# e+e- → η' J/ψ cross section measurement
# arXiv:1911.00885v2 — two η' decay modes
# 11 fb⁻¹ at √s = 4.178–4.600 GeV, 14 energy points

### Dataset description ###
data_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4600"),
]

incMC_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

### Common decay card for e+e- → η' J/ψ, η' → π+π-η (Mode I) and η' → γπ+π- (Mode II) ###
### Both modes share the KKMC psi(4260) top mother ###

# Mode I: η' → π+π-η, η → γγ, J/ψ → ℓ+ℓ-
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 etap J/psi PHSP;
  Enddecay

  Decay etap
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: η' → γπ+π-, J/ψ → ℓ+ℓ-
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 etap J/psi PHSP;
  Enddecay

  Decay etap
  1.0000 gamma pi+ pi- PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC for energy scan ###
exMC_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_etap_jpsi_pipieta"
  config.events = 100_000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_etap_jpsi_gpipi"
  config.events = 100_000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

### ========================================
### Mode I: η' → π+π-η, η → γγ, J/ψ → ℓ+ℓ-
### ========================================
alg_modeI = Algorithm.new("EtapJpsiPipieta")
alg_modeI.set_header(["EtapJpsiPipietaAlg/EtapJpsiPipieta.h"])

sel_modeI = Selection.new
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0
  nGam ">=2"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  identify :pion, against: [:kaon, :proton]
  nlp "==1"
  nlm "==1"
  npip "==1"
  npim "==1"
end
# Reconstruct η → γγ via Kalman fit
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end
# 5C kinematic fit: constrain four-momentum + η mass is already constrained by Kalman
# Participants: ℓ+, ℓ-, π+, π-, η
.kinematic_fit([:lp, :lm, :pip, :pim, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 50
end

alg_modeI.note(:lepton_momentum_cut, "tracks with p > 1 GeV/c assigned as leptons; EMC energy > 1.0 GeV → electron, < 0.4 GeV → muon")
         .note(:jpsi_mass_window, "M(ℓ+ℓ-) required in [3.07, 3.13] GeV/c² (J/ψ signal region); sidebands [3.00,3.06] and [3.14,3.20] used for background")
         .note(:psi2s_veto, "M(η J/ψ) required outside [3.67, 3.70] GeV/c² to veto e+e- → π+π- ψ(2S) → π+π- η J/ψ")
         .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit")
         .note(:photon_conversion_veto, "cosθ(π+π-) < 0.95 to suppress photon conversion background")
         .note(:best_chisq_selection, "if multiple combinations survive, one with smallest χ² is selected")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on(data_points + incMC_points + exMC_modeI)

### ========================================
### Mode II: η' → γπ+π-, J/ψ → ℓ+ℓ-
### ========================================
alg_modeII = Algorithm.new("EtapJpsiGPipi")
alg_modeII.set_header(["EtapJpsiGPipiAlg/EtapJpsiGPipi.h"])

sel_modeII = Selection.new
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0
  nGam ">=1"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  identify :pion, against: [:kaon, :proton]
  nlp "==1"
  nlm "==1"
  npip "==1"
  npim "==1"
end
# 4C kinematic fit for η' → γπ+π- mode
.kinematic_fit([:lp, :lm, :pip, :pim, :gamma]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_modeII.note(:lepton_momentum_cut, "tracks with p > 1 GeV/c assigned as leptons; EMC energy > 1.0 GeV → electron, < 0.4 GeV → muon")
          .note(:jpsi_mass_window, "M(ℓ+ℓ-) required in [3.07, 3.13] GeV/c² (J/ψ signal region); sidebands [3.00,3.06] and [3.14,3.20] used for background")
          .note(:psi2s_veto, "M(π+π- J/ψ) required outside [3.66, 3.71] GeV/c² to veto e+e- → γ_ISR ψ(2S) → γ_ISR π+π- J/ψ")
          .note(:photon_conversion_veto, "cosθ(π+π-) < 0.95 to suppress photon conversion background")
          .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit")
          .note(:best_chisq_selection, "if multiple combinations survive, one with smallest χ² is selected")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on(data_points + incMC_points + exMC_modeII)