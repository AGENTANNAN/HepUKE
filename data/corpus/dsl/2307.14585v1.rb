# BESIII Analysis: Ds+ → μ+ νμ
# Paper: arXiv:2307.14585v1 — Improved measurement of Ds+ → μ+ νμ
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

# Signal MC: e+e- → Ds*± Ds∓, Ds+ → μ+ νμ
# ConExc generator for open-charm processes (continuum + DsDs* production)
decay_card = <<~DECAYCARD
    Decay D_s+
    1.0000  mu+  nu_mu                 ISGW2;
    Enddecay

    Decay D_s-
    1.0000  K+  K-  pi-               PHSP;
    Enddecay

    Decay D_s*+
    1.0000  D_s+  gamma                VSP_PWAVE;
    Enddecay

    Decay D_s*-
    1.0000  D_s-  gamma                VSP_PWAVE;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(ds_data) do |config|
  config.sample_name   = "Ds_munu_signal"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

### Tag analysis ###
alg = TagAnalysis.new("DsMuNu")
alg.set_header(["DsMuNuAlg/DsMuNu.h"])

# ST: Ds- tag with 16 hadronic modes
# 14 modes mapped to authoritative DTagAlg symbols; 2 modes unavailable
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKKPiPi0,
          :DstoPiPiPi,
          :DstoKsK,
          :DstoKsKPi0,
          :DstoKsKsPi,
          :DstoKsKPiPi,
          :DstoEtaGammaGammaPi,
          :DstoEtaPiPiPi0Pi,
          :DstoEtapGammaGammaPiPiPi,
          :DstoEtapGammaPiPiPi,
          :DstoEtaGammaGammaRho,
          :DstoEtaPiPiPi0Rho
  t.charm 1
end

alg.note(:tag_mode_unavailable, "Modes K-π+π- (KPiPi) and ηγγπ+π-π- (EtaGammaGammaPiPiPi) are not available in DTagAlg mode vocabulary; total 16 modes in paper, 14 expressible")

# Signal side: single μ+ + missing νμ
alg.signal_side do |s|
  s.charged(mup: 1, at_least: true)
  s.missing :nu_mu
  s.min_photon_energy 0.025
end

alg.note(:soft_photon, "Transition γ(π0) from Ds* decay selected by minimum |ΔE| method — internal to DTagTool, not DSL-declarable")
alg.note(:muon_pid_emc, "Muon candidate EEMC ∈ (0.0, 0.3) GeV and muon-chamber hit-depth cuts applied in ROOT analysis")
alg.note(:extra_energy_cut, "Emax_extra_γ < 0.3 GeV cut applied in ROOT analysis")

# 4C kinematic fit: four-momentum conservation + Ds/Ds* mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:kinematic_fit_details, "4C fit constrains 4-momentum, Ds masses to known Ds mass, and Ds(*)γ(π0) to nominal Ds* mass; Ds/Ds* mass constraints are internal to tag+fit mechanism")

alg.apply
alg.execute_on(ds_data + ds_incMC + exMC_signal)