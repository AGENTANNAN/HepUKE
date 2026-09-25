# Paper 2308.15362v2: e+e-→K+K-J/ψ cross section at √s 4.61-4.95 GeV
# First observation of Y(4710)
# Hybrid: full reconstruction (both kaons, 4C fit) + partial reconstruction (one kaon missing, 1C fit)
# High-momentum leptons for J/ψ→ℓ+ℓ-
# Multi-energy scan, 12 energy points, 5.85 fb-1

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_points = [
  DatasetManager.real_data.find("706_4610"),  # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),  # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),  # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),  # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),  # 4681.92 MeV
  DatasetManager.real_data.find("706_4700"),  # 4698.82 MeV
  DatasetManager.real_data.find("707_4740"),  # 4739.70 MeV
  DatasetManager.real_data.find("707_4750"),  # 4750.05 MeV
  DatasetManager.real_data.find("707_4780"),  # 4780.54 MeV
  DatasetManager.real_data.find("707_4840"),  # 4843.07 MeV
  DatasetManager.real_data.find("707_4914"),  # 4918.02 MeV
  DatasetManager.real_data.find("707_4946"),  # 4950.93 MeV
]

inc_mc_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# ============================================================
# Full reconstruction channel: e+e- → K+ K- ℓ+ ℓ-
# Both kaons identified, 4C kinematic fit
# ============================================================

decay_card_full = <<~DECAYCARD
  Decay vpho
  1.0000 K+ K- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-  PHSP;
  Enddecay
  End
DECAYCARD

algo_full = Algorithm.new("KKJpsi_full_reco", "00-00-01")
  .set_header(["KKJpsi_full_reco/KKJpsi_full_reco.h"])
  .with_decay_card(decay_card_full)

exMC_full = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KKJpsi_full"
  config.events        = 200_000
  config.decay_card    = decay_card_full
  config.cross_section = :default
end

selection_full = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .pid(method: :probability) do
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==1"
  end
  .kinematic_fit([:kp, :km, :lp, :lm]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_full
  .note(:lepton_separation, "e/μ separation via EMC energy: muon candidates E_EMC < 0.4 GeV, electron candidates E_EMC > 1.0 GeV")
  .note(:bhabha_veto, "ee mode: cosθ(K+K-) < 0.98 to remove radiative Bhabha events")
  .note(:mu_pi_misid, "μμ mode: at least one muon MUC penetration depth > 30 cm to suppress μ/π misidentification")
  .apply(selection_full)

# ============================================================
# Partial reconstruction channel: e+e- → K± K∓_miss ℓ+ ℓ-
# One kaon missing, 1C kinematic fit
# ============================================================

decay_card_partial = <<~DECAYCARD
  Decay vpho
  1.0000 K+ K- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-  PHSP;
  Enddecay
  End
DECAYCARD

algo_partial = Algorithm.new("KKJpsi_partial_reco", "00-00-01")
  .set_header(["KKJpsi_partial_reco/KKJpsi_partial_reco.h"])
  .with_decay_card(decay_card_partial)

exMC_partial = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KKJpsi_partial"
  config.events        = 200_000
  config.decay_card    = decay_card_partial
  config.cross_section = :default
end

selection_partial = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .pid(method: :probability) do
    identify :kaon, against: [:pion]
  end
  .kinematic_fit([:lp, :lm, :kp]) do
    miss_track_of :km
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_partial
  .note(:lepton_separation, "e/μ separation via EMC energy: muon candidates E_EMC < 0.4 GeV, electron candidates E_EMC > 1.0 GeV")
  .note(:bhabha_veto_partial, "ee mode: cosθ(e+) < 0.8, cosθ(e-) > -0.8, |cosθ(K±)| < 0.8; also require |cos(α_K±K_miss)| < 0.95, |cos(α_K±e±)| < 0.95, |cos(α_K±e∓)| < 0.95 to suppress radiative Bhabha with γ conversion")
  .note(:mu_pi_misid_partial, "μμ mode: both muon candidates required to have MUC penetration depth > 30 cm")
  .note(:momentum_cut, "lepton p > 0.95 GeV/c for √s < 4.84 GeV, p > 1.05 GeV/c for √s ≥ 4.84 GeV")
  .apply(selection_partial)

# Execute on datasets
algo_full.execute_on(data_points + inc_mc_points + exMC_full)
algo_partial.execute_on(data_points + inc_mc_points + exMC_partial)