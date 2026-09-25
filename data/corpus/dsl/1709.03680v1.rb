# Paper: 1709.03680v1 — Measurements of branching fractions for semi-leptonic decays
#   Ds+ → φ e+ νe, φ μ+ νμ, η μ+ νμ, η' μ+ νμ
# Analysis type: TagAnalysis — ST/DT with semileptonic signal side
# Data: √s = 4.009 GeV, 482 pb⁻¹

# Load dataset tables
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# Data: √s = 4.009 GeV
# ============================================================
data_4009 = DatasetManager.real_data.find("703_4009")
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")

# ============================================================
# Decay card for ψ(4040) → DsDs at 4.009 GeV
# ============================================================
decay_card_dsds = <<~DECAYCARD
  Decay psi(4040)
  1.0 D_s+ D_s- VSS;
  Enddecay
  End
DECAYCARD

# ============================================================
# Common tag side: Ds- hadronic decay (10 modes)
# ============================================================

# ============================================================
# Mode 1: Ds+ → φ e+ νe
# ============================================================
alg_phi_e = TagAnalysis.new("DsPhiENu")

alg_phi_e.set_header(["DsPhiENuAlg/DsPhiENu.h"])
          .set_constant({ "ECMS" => [:double, 4.009] })
          .with_decay_card(decay_card_dsds)

alg_phi_e.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg_phi_e.signal_side do |s|
  s.charged(kp: 1, km: 1, ep: 1)
  s.require_charge 0           # K+ K- e+ = +1+(-1)+1 = 1, but Ds+=+1... wait
  # Ds+ → φ e+ νe: φ→K+K- gives 0 net charge from φ, +1 from e+, total +1
  s.require_charge 1
  s.missing :nu_e              # massless neutrino
end

alg_phi_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :km).constrain_to_nominal_mass_of(:phi)
  f.chi2_cut 200
end

alg_phi_e.apply

# ============================================================
# Mode 2: Ds+ → φ μ+ νμ
# ============================================================
alg_phi_mu = TagAnalysis.new("DsPhiMuNu")

alg_phi_mu.set_header(["DsPhiMuNuAlg/DsPhiMuNu.h"])
           .set_constant({ "ECMS" => [:double, 4.009] })
           .with_decay_card(decay_card_dsds)

alg_phi_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg_phi_mu.signal_side do |s|
  s.charged(kp: 1, km: 1, mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_phi_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :km).constrain_to_nominal_mass_of(:phi)
  f.chi2_cut 200
end

alg_phi_mu.apply

# ============================================================
# Mode 3: Ds+ → η μ+ νμ  (η → γγ)
# ============================================================
alg_eta_mu = TagAnalysis.new("DsEtaMuNu")

alg_eta_mu.set_header(["DsEtaMuNuAlg/DsEtaMuNu.h"])
           .set_constant({ "ECMS" => [:double, 4.009] })
           .with_decay_card(decay_card_dsds)

alg_eta_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg_eta_mu.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_eta_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_eta_mu.apply

# ============================================================
# Mode 4: Ds+ → η' μ+ νμ  (η' → ηπ+π-, η → γγ)
# ============================================================
alg_etap_mu = TagAnalysis.new("DsEtapMuNu")

alg_etap_mu.set_header(["DsEtapMuNuAlg/DsEtapMuNu.h"])
            .set_constant({ "ECMS" => [:double, 4.009] })
            .with_decay_card(decay_card_dsds)

alg_etap_mu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg_etap_mu.signal_side do |s|
  s.photons 2                    # for η → γγ
  s.charged(pip: 1, pim: 1, mup: 1)  # π+π- from η'→ηπ+π-, μ+ from SL
  s.require_charge 1
  s.missing :nu_mu
end

alg_etap_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  # η' mass constraint: ηπ+π- system; the eta is already constrained
  # DSL limitation: intermediate constraint nesting not supported in v1
  f.chi2_cut 200
end

alg_etap_mu.apply

# ============================================================
# Execute — note: four separate TagAnalysis objects
# ============================================================
sig_mc_ds = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_ds_semilep"
  config.related_dataset = data_4009
  config.events          = 500_000
  config.decay_card      = decay_card_dsds
  config.cross_section   = :default
end

alg_phi_e.execute_on([data_4009, incMC_4009, sig_mc_ds])
alg_phi_mu.execute_on([data_4009, incMC_4009, sig_mc_ds])
alg_eta_mu.execute_on([data_4009, incMC_4009, sig_mc_ds])
alg_etap_mu.execute_on([data_4009, incMC_4009, sig_mc_ds])