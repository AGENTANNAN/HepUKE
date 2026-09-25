# DSL: e+e- → γ η(') ηc  at √s = 4.258–4.681 GeV (multi-energy R-scan)
# Paper: 2504.13539v1 — Search for 1-+ charmonium-like hybrid
# ηc reconstructed via 16 hadronic modes (representative mode: p pbar used here)
# η' reconstructed via η π+ π- , η → γγ
# NOTE: 16 ηc decay modes handled simultaneously at ROOT level — DSL expresses
# one representative mode per process; other modes require analogous fit declarations.

# ============================================================
# Datasets — representative energy point √s = 4.681 GeV
# ============================================================
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4681  = DatasetManager.real_data.find("708_4681")
incMC_4681 = DatasetManager.inclusive_mc.find("708_4681")

# ============================================================
# Decay card — KKMC + psi(4260) top mother (continuum production)
# Signal: e+e- → γ η ηc  and  e+e- → γ η' ηc
# ============================================================
decay_card_gam_eta_etac = <<~DECAYCARD
  Decay psi(4260)
  1 gamma eta eta_c PHSP;
  Enddecay
  Decay eta_c
  1 p+ anti-p- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

decay_card_gam_etap_etac = <<~DECAYCARD
  Decay psi(4260)
  1 gamma eta' eta_c PHSP;
  Enddecay
  Decay eta_c
  1 p+ anti-p- PHSP;
  Enddecay
  Decay eta'
  1 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Signal MC — single energy point (representative)
# ============================================================
sig_gam_eta_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_gam_eta_etac_4681"
  config.related_dataset = data_4681
  config.events          = 100_000
  config.decay_card      = decay_card_gam_eta_etac
  config.cross_section   = :default
end

sig_gam_etap_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_gam_etap_etac_4681"
  config.related_dataset = data_4681
  config.events          = 100_000
  config.decay_card      = decay_card_gam_etap_etac
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: e+e- → γ η ηc  (ηc → p pbar)
# Final state: γ (radiative) + η (→γγ) + p + pbar
# ============================================================
alg1 = Algorithm.new("GamEtaEtac")
alg1.set_header(["GamEtaEtacAlg/GamEtaEtac.h"])
    .set_constant({ "ECMS" => [:double, 4.681] })

sel1 = Selection.new
  # Charged tracks: p, pbar
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
  end
  # Photon selection: at least 3 photons (radiative + η→γγ)
  .select_photon do
    nGam ">=3"
    min_energy 0.025
    min_angle 10.0
  end
  # PID — probability method
  .pid do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct η → γγ via Kalman fit (1C mass constraint)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing-hypothesis veto: 2γ + hadrons (no η→γγ)
  .kinematic_fit([:gamma, :gamma, :prp, :prm]) do
    constrain_four_momentum
  end
  # Competing-hypothesis veto: 4γ + hadrons
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prp, :prm]) do
    constrain_four_momentum
  end
  # Nominal 4C fit: radiative γ + η + p + pbar
  .kinematic_fit([:gamma, :eta, :prp, :prm]) do
    constrain_four_momentum
    chi2_cut 20
    nominal
  end

alg1.with_decay_card(decay_card_gam_eta_etac).apply(sel1)
alg1.note(:multi_mode, "16 eta_c hadronic modes: p pbar, 2(pi+pi-), 2(K+K-), K+K-pi+pi-, p pbar pi+pi-, 3(pi+pi-), K+K-2(pi+pi-), K+K-pi0, p pbar pi0, K_S0 K± pi∓, K_S0 K± pi∓ pi± pi∓, pi+pi-eta, K+K-eta, 2(pi+pi-)eta, pi+pi-pi0pi0, 2(pi+pi-)pi0pi0. ROOT handles mode combination and best-candidate selection via min(χ²_tot).")
alg1.note(:competing_hypothesis, "χ²_4C(3γ+hadrons) < χ²_4C(2γ+hadrons) and χ²_4C(4γ+hadrons) — evaluated in ROOT by comparing stored χ² values.")
alg1.note(:ks_reconstruction, "For eta_c modes containing K_S0: reconstruct via secondary vertex fit (pip pim) with mass window |M_ππ - m_K_S0| < 20 MeV and decay length > 2σ.")
alg1.note(:pi0_reconstruction, "For eta_c modes containing π0: reconstruct via Kalman fit (γγ → π0) with mass window |M_γγ - m_π0| < 15 MeV.")
alg1.execute_on([data_4681, incMC_4681, sig_gam_eta_etac])

# ============================================================
# Algorithm 2: e+e- → γ η' ηc  (η' → η π+ π- , ηc → p pbar)
# Final state: γ (radiative) + η (from η') + π+ + π- + p + pbar
# ============================================================
alg2 = Algorithm.new("GamEtapEtac")
alg2.set_header(["GamEtapEtacAlg/GamEtapEtac.h"])
    .set_constant({ "ECMS" => [:double, 4.681] })

sel2 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
  end
  .select_photon do
    nGam ">=3"
    min_energy 0.025
    min_angle 10.0
  end
  .pid do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct η → γγ (from η' decay)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing-hypothesis veto: 2γ + π+π- + hadrons
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :prp, :prm]) do
    constrain_four_momentum
  end
  # Competing-hypothesis veto: 4γ + π+π- + hadrons
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :prp, :prm]) do
    constrain_four_momentum
  end
  # Nominal 4C fit: radiative γ + η + π+ + π- + p + pbar
  # χ²_4C(3γ π+π- + hadrons) < 30
  .kinematic_fit([:gamma, :eta, :pip, :pim, :prp, :prm]) do
    constrain_four_momentum
    chi2_cut 30
    nominal
  end

alg2.with_decay_card(decay_card_gam_etap_etac).apply(sel2)
alg2.note(:multi_mode, "Same 16 eta_c hadronic modes as Algorithm 1. η' reconstructed via η π+ π- . ROOT handles combined mass-difference criterion for η_c/η' candidate selection.")
alg2.note(:competing_hypothesis, "χ²_4C(3γ π+π- + hadrons) < χ²_4C(2γ π+π- + hadrons) and χ²_4C(4γ π+π- + hadrons) — evaluated in ROOT.")
alg2.note(:ks_reconstruction, "For eta_c modes containing K_S0: secondary vertex fit as in Algorithm 1.")
alg2.note(:pi0_reconstruction, "For eta_c modes containing π0: Kalman fit (γγ → π0) as in Algorithm 1.")
alg2.execute_on([data_4681, incMC_4681, sig_gam_etap_etac])