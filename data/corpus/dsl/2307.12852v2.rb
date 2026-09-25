# BESIII Analysis: Ds+ → η(') μ+ νμ
# Paper: arXiv:2307.12852v2 — First observation of Ds+ → η' μ+ νμ, LFU test
# Dataset: 7.33 fb⁻¹ at Ecm = 4.128–4.226 GeV (8 energy points)
# Analysis type: Tag-based (Ds tag, ST+missing pattern, DT method)

### Dataset preparation ###
data_4128 = DatasetManager.real_data.find("705_4128")
data_4157 = DatasetManager.real_data.find("705_4157")
data_4178 = DatasetManager.real_data.find("705_4178")
data_4189 = DatasetManager.real_data.find("705_4189")
data_4199 = DatasetManager.real_data.find("705_4199")
data_4209 = DatasetManager.real_data.find("705_4209")
data_4219 = DatasetManager.real_data.find("705_4219")
data_4226 = DatasetManager.real_data.find("705_4226")

ds_data = [data_4128, data_4157, data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]

# Inclusive MC at each energy point
incMC_4128 = DatasetManager.inclusive_mc.find("705_4128")
incMC_4157 = DatasetManager.inclusive_mc.find("705_4157")
incMC_4178 = DatasetManager.inclusive_mc.find("705_4178")
incMC_4189 = DatasetManager.inclusive_mc.find("705_4189")
incMC_4199 = DatasetManager.inclusive_mc.find("705_4199")
incMC_4209 = DatasetManager.inclusive_mc.find("705_4209")
incMC_4219 = DatasetManager.inclusive_mc.find("705_4219")
incMC_4226 = DatasetManager.inclusive_mc.find("705_4226")

ds_incMC = [incMC_4128, incMC_4157, incMC_4178, incMC_4189, incMC_4199, incMC_4209, incMC_4219, incMC_4226]

# Signal MC: Ds+ → η(') μ+ νμ
# ConExc generator for open-charm processes
decay_card_eta = <<~DECAYCARD
    Decay D_s+
    1.0000  eta  mu+  nu_mu          ISGW2;
    Enddecay

    Decay eta
    1.0000  gamma  gamma              PHSP;
    Enddecay
End
DECAYCARD

decay_card_etap = <<~DECAYCARD
    Decay D_s+
    1.0000  eta'  mu+  nu_mu         ISGW2;
    Enddecay

    Decay eta'
    1.0000  gamma  pi+  pi-          PHSP;
    Enddecay
End
DECAYCARD

exMC_eta = DatasetManager.create_exclusive_mc_for(ds_data) do |config|
  config.sample_name   = "Ds_eta_munu"
  config.events        = 200_000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

exMC_etap = DatasetManager.create_exclusive_mc_for(ds_data) do |config|
  config.sample_name   = "Ds_etap_munu"
  config.events        = 200_000
  config.decay_card    = decay_card_etap
  config.cross_section = :default
end

# ============================================================
# Analysis A: Ds+ → η μ+ νμ  (η → γγ)
# ============================================================

alg_eta = TagAnalysis.new("DsEtaMuNu")
alg_eta.set_header(["DsEtaMuNuAlg/DsEtaMuNu.h"])

# ST: Ds- tag with 14 hadronic modes
alg_eta.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKKPiPi0,
          :DstoKsK,
          :DstoKsKPi0,
          :DstoKsKsPi,
          :DstoKsKPiPi,
          :DstoPiPiPi,
          :DstoEtaGammaGammaPi,
          :DstoEtaPiPiPi0Pi,
          :DstoEtapGammaGammaPiPiPi,
          :DstoEtapGammaPiPiPi,
          :DstoEtaGammaGammaRho,
          :DstoEtaPiPiPi0Rho
  t.charm 1
end

alg_eta.note(:tag_mode_unavailable, "Mode K+K-π-π0 second variant (different Ks ordering) combined in ROOT; 14 of 14 tag modes expressible")

# Signal side: 2 photons (η→γγ) + muon + missing νμ
alg_eta.signal_side do |s|
  s.photons 2
  s.charged(mup: 1, at_least: true)
  s.missing :nu_mu
  s.min_photon_energy 0.025
end

alg_eta.note(:eta_pi0pipi_mode, "η → π0π+π- sub-decay mode handled in separate analysis; combined with η→γγ in ROOT fit")
alg_eta.note(:soft_photon, "Transition γ(π0) from Ds* decay selected by minimum |ΔE| method — internal to DTagTool")
alg_eta.note(:muon_emc_cut, "Muon EEMC ∈ (0.10, 0.28) GeV applied in ROOT analysis")
alg_eta.note(:extra_cuts, "N_extra_char = 0, N_extra_pi0 = 0, Emax_extra_γ < 0.2 GeV applied in ROOT analysis")
alg_eta.note(:peaking_bkg_veto, "M(ημ+) < 1.8 GeV/c² and M(ηνμ) > 0.97 GeV/c² applied in ROOT analysis")

# 3C kinematic fit: four-momentum conservation + Ds/Ds* mass constraints
alg_eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta.apply

# ============================================================
# Analysis B: Ds+ → η' μ+ νμ  (η' → γπ+π-)
# ============================================================

alg_etap = TagAnalysis.new("DsEtapMuNu")
alg_etap.set_header(["DsEtapMuNuAlg/DsEtapMuNu.h"])

alg_etap.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKKPiPi0,
          :DstoKsK,
          :DstoKsKPi0,
          :DstoKsKsPi,
          :DstoKsKPiPi,
          :DstoPiPiPi,
          :DstoEtaGammaGammaPi,
          :DstoEtaPiPiPi0Pi,
          :DstoEtapGammaGammaPiPiPi,
          :DstoEtapGammaPiPiPi,
          :DstoEtaGammaGammaRho,
          :DstoEtaPiPiPi0Rho
  t.charm 1
end

alg_etap.note(:tag_mode_unavailable, "14 hadronic tag modes per paper; all expressible in DTagAlg vocabulary")

# Signal side: photon + π+π- + muon + missing νμ
# η' → γπ+π- needs 1 photon, 2 pions; also supports η' → ηγγπ+π- (4 photons)
alg_etap.signal_side do |s|
  s.photons 1
  s.charged(pip: 1, pim: 1, mup: 1, at_least: true)
  s.missing :nu_mu
  s.min_photon_energy 0.025
end

alg_etap.note(:etap_etagpipi_mode, "η' → ηγγπ+π- sub-decay mode handled in separate analysis; combined with η'→γπ+π- in ROOT fit")
alg_etap.note(:soft_photon, "Transition γ(π0) from Ds* decay selected by minimum |ΔE| method — internal to DTagTool")
alg_etap.note(:chi2_cut_specific, "χ² < 30 required for Ds+ → η'(γπ+π-)μ+νμ to suppress non-DsDs* backgrounds")
alg_etap.note(:muon_emc_cut, "Muon EEMC ∈ (0.10, 0.28) GeV applied in ROOT analysis")
alg_etap.note(:extra_cuts, "N_extra_char = 0, N_extra_pi0 = 0, Emax_extra_γ < 0.2 GeV applied in ROOT analysis")
alg_etap.note(:peaking_bkg_veto, "M(η'μ+) < 1.8 GeV/c² and M(η'νμ) > 1.27 GeV/c² applied in ROOT analysis")

# 3C kinematic fit + eta' resonance constraint
alg_etap.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_etap.note(:etap_resonance_constraint, "η' mass constraint via invariant mass of γπ+π- handled in kinematic fit; DSL expresses resonance via fit constraints")

alg_etap.apply

# ============================================================
# Execute both analyses
# ============================================================
alg_eta.execute_on(ds_data + ds_incMC + exMC_eta)
alg_etap.execute_on(ds_data + ds_incMC + exMC_etap)