# Paper: 2203.05815v2 — Observation of resonance structures in e+e- → π+π- ψ2(3823)
# and mass measurement of ψ2(3823)
# Multi-energy scan: 4.23–4.70 GeV, ORDINARY analysis (cross-section scan)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ── Datasets: 20 energy-scan points ──────────────────────────────────
data_points = [
  DatasetManager.real_data.find("703_4230"),   # 4.2263 GeV
  DatasetManager.real_data.find("703_4260"),   # 4.2580 GeV
  DatasetManager.real_data.find("705_4290"),   # 4.2879 GeV
  DatasetManager.real_data.find("705_4315"),   # 4.3121 GeV
  DatasetManager.real_data.find("705_4340"),   # 4.3374 GeV
  DatasetManager.real_data.find("703_4360"),   # 4.3583 GeV
  DatasetManager.real_data.find("705_4380"),   # 4.3774 GeV
  DatasetManager.real_data.find("705_4400"),   # 4.3965 GeV
  DatasetManager.real_data.find("703_4420"),   # 4.4156 GeV
  DatasetManager.real_data.find("705_4440"),   # 4.4362 GeV
  DatasetManager.real_data.find("703_4470"),   # 4.4671 GeV
  DatasetManager.real_data.find("703_4530"),   # 4.5271 GeV
  DatasetManager.real_data.find("703_4575"),   # 4.5745 GeV
  DatasetManager.real_data.find("703_4600"),   # 4.5995 GeV
  DatasetManager.real_data.find("706_4610"),   # 4.6120 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.6278 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.6408 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.6613 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.6811 GeV
  DatasetManager.real_data.find("706_4700"),   # 4.6984 GeV
]

# Inclusive MC for each energy point
inc_mc_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# ── Decay card ────────────────────────────────────────────────────────
# e+e- → π+π- ψ2(3823), ψ2(3823) → γ χc1, χc1 → γ J/ψ, J/ψ → ℓ+ℓ−
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3823)  PHSP;
  Enddecay
  Decay psi(3823)
  1.0000 gamma chi_c1  PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.0000 l+ l-  PHSP;
  Enddecay
  End
DECAYCARD

# ── Exclusive MC ──────────────────────────────────────────────────────
ex_mc_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pipi_psi23823"
  config.events        = 50000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# ── Algorithm ─────────────────────────────────────────────────────────
algorithm = Algorithm.new("PipiPsi23823Scan")

# ── Event selection ───────────────────────────────────────────────────
# 4 charged tracks, net charge zero; at least 1 photon (partial reco) or ≥2 (full reco)
# Lepton PID for J/ψ → ℓ+ℓ−; both e+e- and μ+μ- channels combined
event_selection = Selection.new
  .select_track do
    nChrp "==2"
    nChrn "==2"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=1"
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:ep, :em, against: :pion)
    identify(:mup, :mum, against: :pion)
    prob_cut 0.001
  end
  .assign(lp: :ep, lm: :em)
  .kinematic_fit([:pip, :pim, :lp, :lm, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algorithm
  .note(:lepton_universality, "both e+e- and μ+μ- channels are used for J/ψ reconstruction; lepton candidate lists from both electron and muon PID are combined; the invariant mass window 3.06 < M(ℓ+ℓ-) < 3.135 GeV/c2 selects J/ψ candidates with sideband subtraction in ROOT")
  .note(:helix_correction, "track helix parameter correction applied to each MC-simulated event during the 1C kinematic fit (Ref.[40] M.Ablikim et al., PRD 87, 012002); 1.7% systematic uncertainty assigned")
  .note(:partial_reconstruction, "two reconstruction strategies: (i) partial reco with Nγ=1 and a missing photon inferred from recoil, 1C kinematic fit; (ii) full reco with Nγ≥2 and 4C kinematic fit; the signal extraction uses a simultaneous fit to M_recoil(π+π-) in both categories in ROOT")
  .note(:background_veto, "cos(opening_angle(π+,π-)) < 0.98 (gamma-conversion veto); M(γγ_miss π+π-) > 0.65 GeV/c2 (ηJ/ψ veto); |M(π+π-J/ψ)−m[ψ(2S)]| > 7 MeV/c2 (ψ(2S) veto using M(π+π-ℓ+ℓ-)−M(ℓ+ℓ-)+m(J/ψ))")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on(data_points + inc_mc_points + ex_mc_signal)