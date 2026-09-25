### Datasets: ψ(4040) at √s = 4.009 GeV (482 pb⁻¹) ###
data_4009  = DatasetManager.real_data.find("703_4009")     # 4.009 GeV real data (482 pb⁻¹)
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")  # Matching inclusive MC

# Decay card for the signal process (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(4040)
    1.0000 D_s+ D_s- VSS;
    Enddecay
    End
DECAYCARD

# Exclusive MC: ψ(4040) → Ds+ Ds-, 500k events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4009_DsDs"
  config.related_dataset = data_4009
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (single-tag Ds- → K+ K- π-, signal side semileptonic) ###

# ------------------------------------------------------------------
# Mode 1: Ds+ → φ e+ νe  (φ → K+ K-)
# ------------------------------------------------------------------
alg_phi_e = TagAnalysis.new("DsToPhiENu")
alg_phi_e.set_header(["DsToPhiENuAlg/DsToPhiENu.h"])
         .set_constant({"ECMS" => [:double, 4.009]})
         .with_decay_card(decay_card_signal)

alg_phi_e.tag_side(:Ds) do |t|
  t.modes :DstoKKPi   # tag Ds- → K+ K- π-
  t.charm -1          # pin the tagged side to Ds-
end

alg_phi_e.signal_side do |s|
  s.charged(kp: 1, km: 1, ep: 1)  # K+ K- e+ from φ e+ νe
  s.require_charge 1              # net charge +1
  s.missing :nu_e                 # one missing massless νe
end

alg_phi_e.fit do |f|
  f.constrain_four_momentum                                        # 4-momentum conservation
  f.invariant_mass_of(:kp, :km).constrain_to_nominal_mass_of(:phi) # φ mass constraint on K+K-
  f.chi2_cut 200
end

alg_phi_e.apply                       # tag spec — takes no Selection argument
alg_phi_e.execute_on([data_4009, incMC_4009, exMC_signal])

# ------------------------------------------------------------------
# Mode 2: Ds+ → φ μ+ νμ  (φ → K+ K-)
# ------------------------------------------------------------------
alg_phi_mu = TagAnalysis.new("DsToPhiMuNu")
alg_phi_mu.set_header(["DsToPhiMuNuAlg/DsToPhiMuNu.h"])
          .set_constant({"ECMS" => [:double, 4.009]})
          .with_decay_card(decay_card_signal)

alg_phi_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi   # tag Ds- → K+ K- π-
  t.charm -1
end

alg_phi_mu.signal_side do |s|
  s.charged(kp: 1, km: 1, mup: 1)  # K+ K- μ+ from φ μ+ νμ
  s.require_charge 1               # net charge +1
  s.missing :nu_mu                 # one missing massless νμ
end

alg_phi_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :km).constrain_to_nominal_mass_of(:phi) # φ mass constraint on K+K-
  f.chi2_cut 200
end

alg_phi_mu.apply
alg_phi_mu.execute_on([data_4009, incMC_4009, exMC_signal])

# ------------------------------------------------------------------
# Mode 3: Ds+ → η μ+ νμ  (η → γγ)
# ------------------------------------------------------------------
alg_eta_mu = TagAnalysis.new("DsToEtaMuNu")
alg_eta_mu.set_header(["DsToEtaMuNuAlg/DsToEtaMuNu.h"])
          .set_constant({"ECMS" => [:double, 4.009]})
          .with_decay_card(decay_card_signal)

alg_eta_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi   # tag Ds- → K+ K- π-
  t.charm -1
end

alg_eta_mu.signal_side do |s|
  s.photons 2             # two photons from η → γγ
  s.charged(mup: 1)       # μ+ from η μ+ νμ
  s.require_charge 1      # net charge +1
  s.missing :nu_mu        # one missing massless νμ
end

alg_eta_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta) # η mass constraint on γγ
  f.chi2_cut 200
end

alg_eta_mu.apply
alg_eta_mu.execute_on([data_4009, incMC_4009, exMC_signal])

# ------------------------------------------------------------------
# Mode 4: Ds+ → η' μ+ νμ  (η' → η π+ π-, η → γγ)
# ------------------------------------------------------------------
alg_etap_mu = TagAnalysis.new("DsToEtaPrimeMuNu")
alg_etap_mu.set_header(["DsToEtaPrimeMuNuAlg/DsToEtaPrimeMuNu.h"])
           .set_constant({"ECMS" => [:double, 4.009]})
           .with_decay_card(decay_card_signal)

alg_etap_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi   # tag Ds- → K+ K- π-
  t.charm -1
end

alg_etap_mu.signal_side do |s|
  s.photons 2                        # two photons from η → γγ
  s.charged(pip: 1, pim: 1, mup: 1)  # π+ π- from η'→ηπ+π-, μ+ from η' μ+ νμ
  s.require_charge 1                 # net charge +1
  s.missing :nu_mu                   # one missing massless νμ
end

alg_etap_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta) # η mass constraint on γγ
  f.chi2_cut 200
end

alg_etap_mu.apply
alg_etap_mu.execute_on([data_4009, incMC_4009, exMC_signal])