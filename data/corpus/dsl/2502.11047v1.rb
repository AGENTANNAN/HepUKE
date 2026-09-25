# BESIII DSL: Search for Lambda_c+ -> Sigma0 K+ pi0 and Sigma0 K+ pi+ pi-
# Paper: 2502.11047v1
# CMS energy: 4599.53-4698.82 MeV (7 energy points, 4.5 fb^-1)
# Method: Single-tag (reconstruct one Lambda_c+, no recoil requirement)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 7 energy points for Lambda_c+ Lambda_c- production above threshold
ALL_ENERGY_KEYS = %w[703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700]

all_data = ALL_ENERGY_KEYS.map { |k| DatasetManager.real_data.find(k) }
all_inc_mc = ALL_ENERGY_KEYS.map { |k| DatasetManager.inclusive_mc.find(k) }

# ============================================================
# Decay cards
# Lambda_c+ -> Sigma0 K+ pi0, Sigma0 -> Lambda gamma, Lambda -> p pi-, pi0 -> gamma gamma
# ============================================================

decay_card_sig0kpi0 = <<~DECAYCARD
  Decay Lambda_c+
  1.000 Sigma0 K+ pi0 PHSP;
  Enddecay
  Decay Sigma0
  1.000 Lambda gamma PHSP;
  Enddecay
  Decay Lambda
  1.000 p pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

decay_card_sig0kpipipi = <<~DECAYCARD
  Decay Lambda_c+
  1.000 Sigma0 K+ pi+ pi- PHSP;
  Enddecay
  Decay Sigma0
  1.000 Lambda gamma PHSP;
  Enddecay
  Decay Lambda
  1.000 p pi- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive signal MC (7 energy points for each mode)
sig_mc_pi0 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_Lc_Sig0Kpi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card_sig0kpi0
  config.cross_section = :default
end

sig_mc_pipipi = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_Lc_Sig0Kpipipi"
  config.events        = 1_000_000
  config.decay_card    = decay_card_sig0kpipipi
  config.cross_section = :default
end

# ============================================================
# Algorithm 1: Lambda_c+ -> Sigma0 K+ pi0
# ============================================================

alg_pi0 = Algorithm.new("LcToSig0Kpi0")
alg_pi0.set_header(["LcToSig0Kpi0/LcToSig0Kpi0.h"])
alg_pi0.set_constant({ "ECMS" => [:double, 4.600] })
alg_pi0.with_decay_card(decay_card_sig0kpi0)

evt_sel_pi0 = Selection.new
  # Charged tracks: at least one K+, one pi0 from Lambda decay (p, pi-)
  .select_track do
    nChrp ">=1"   # K+ and proton from Lambda
    nChrn ">=1"   # pi- from Lambda
    nChr ">=2"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  # PID: proton (L(p) > L(K) and L(p) > L(pi)), kaon (L(K) > L(pi)), pion (L(pi) > L(K))
  .pid do
    identify :prp, :kp, :pim
  end
  # Photons: at least 3 (2 for pi0, 1 for Sigma0 -> Lambda gamma)
  .select_photon do
    nGam ">=3"
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  end
  # Photon timing (EMC time within [0, 700] ns) and angle > 10 deg from charged tracks
  # Lambda reconstruction: p pi- with vertex fit
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda)
    chi2_cut 100
    require_decay_length 2.0  # > 2 * vertex resolution
  end
  # Lambda mass window: [1.111, 1.121] GeV/c^2
  .where("abs(M(p pi-) - 1.1157) < 0.005")

# pi0 reconstruction from gamma gamma
evt_sel_pi0
  .where("M(gamma gamma) > 0.115 && M(gamma gamma) < 0.150") do
    # pi0 mass window before 1C fit
  end
  # 1C kinematic fit constraining gamma gamma to pi0 mass
  .kinematic_fit([:gamma, :gamma]) do
    constrain_invariant_mass_of([:gamma, :gamma]).to_nominal(:pi0)
    chi2_cut 200
    nominal
  end

# Sigma0 reconstruction from Lambda + gamma
evt_sel_pi0
  .where("M(Lambda gamma) > 1.179 && M(Lambda gamma) < 1.203") do
    # Sigma0 mass window [1.179, 1.203] GeV
  end

# Opening angle between photon and antiproton > 20 deg (suppress p EMC noise)
evt_sel_pi0
  .where("cos(angle(p_bar, gamma)) < cos(20*pi/180)") do
    # theta_{p_bar gamma} > 20 deg
  end

# DeltaE suppression cut for peaking background (Xi0 K+, Lambda K*+)
evt_sel_pi0
  .where("DeltaE_p_pi-_K+_gamma_gamma > -0.160 && DeltaE_p_pi-_K+_gamma_gamma < -0.030") do
    # -160 < DeltaE(p pi- K+ gamma gamma) < -30 MeV
  end

# Lambda_c+ candidate: best combination by minimum |DeltaE|
# DeltaE window: [-27, 6] MeV
# M_BC signal region: [2.282, 2.291] GeV/c^2

alg_pi0.apply(evt_sel_pi0)

alg_pi0.note(:single_tag, "Single-tag method: reconstruct one Lambda_c+ with no recoil requirement")
alg_pi0.note(:lambda_veto, "Lambda candidates with |Vz| < 20 cm for daughter proton")
alg_pi0.note(:pi0_1c_fit, "1C kinematic fit constraining M(gamma gamma) to pi0 mass, chi2 < 200")
alg_pi0.note(:sigma0_mass, "Sigma0 mass window: [1.179, 1.203] GeV/c^2")
alg_pi0.note(:deltaE_cut, "DeltaE window: [-27, 6] MeV applied in ROOT from stored variable")
alg_pi0.note(:mbc_signal, "M_BC signal region: [2.282, 2.291] GeV/c^2 applied in ROOT")
alg_pi0.note(:best_combination, "Best Lambda_c+ candidate chosen by minimum |DeltaE|")
alg_pi0.note(:peaking_bkg_suppression, "DeltaE(p pi- K+ gamma gamma) in [-160, -30] MeV suppresses Xi0 K+ and Lambda K*+ bkg")
alg_pi0.note(:theta_pbar_gamma, "Opening angle between anti-proton and photon > 20 deg")
alg_pi0.note(:multi_energy, "7 energy points 4599.53-4698.82 MeV fitted simultaneously; ARGUS bkg + signal MC convolved with Gaussian")
alg_pi0.note(:upper_limit, "Upper limit at 90% CL determined via likelihood scan with systematics convolved")

alg_pi0.execute_on(all_data + all_inc_mc + [sig_mc_pi0])

# ============================================================
# Algorithm 2: Lambda_c+ -> Sigma0 K+ pi+ pi-
# ============================================================

alg_pipipi = Algorithm.new("LcToSig0Kpipipi")
alg_pipipi.set_header(["LcToSig0Kpipipi/LcToSig0Kpipipi.h"])
alg_pipipi.set_constant({ "ECMS" => [:double, 4.600] })
alg_pipipi.with_decay_card(decay_card_sig0kpipipi)

evt_sel_pipipi = Selection.new
  # Charged tracks: K+, pi+, pi- from signal, p, pi- from Lambda
  .select_track do
    nChrp ">=2"   # K+ and proton
    nChrn ">=2"   # 2 pi- (one from Lambda, one from signal)
    nChr ">=4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  # PID: proton, kaon, pions
  .pid do
    identify :prp, :kp, :pip, :pim
  end
  # Photons: at least 1 (for Sigma0 -> Lambda gamma; no pi0 in this mode)
  .select_photon do
    nGam ">=1"
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  end
  # Lambda reconstruction
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda)
    chi2_cut 100
    require_decay_length 2.0
  end
  .where("abs(M(p pi-) - 1.1157) < 0.005")

# Sigma0 reconstruction from Lambda + gamma
evt_sel_pipipi
  .where("M(Lambda gamma) > 1.179 && M(Lambda gamma) < 1.203") do
    # Sigma0 mass window
  end

# DeltaE suppression: DeltaE(p pi- K+ pi+ pi-) < -40 MeV
evt_sel_pipipi
  .where("DeltaE_p_pi-_K+_pi+_pi- < -0.040") do
    # peaking background suppression
  end

# Lambda_c+ candidate: best combination by minimum |DeltaE|
# DeltaE window: [-21, 7] MeV
# M_BC signal region: [2.282, 2.291] GeV/c^2

alg_pipipi.apply(evt_sel_pipipi)

alg_pipipi.note(:single_tag, "Single-tag method: reconstruct one Lambda_c+ with no recoil requirement")
alg_pipipi.note(:lambda_reco, "Lambda -> p pi- with secondary vertex fit chi2 < 100")
alg_pipipi.note(:sigma0_mass, "Sigma0 mass window: [1.179, 1.203] GeV/c^2")
alg_pipipi.note(:deltaE_cut, "DeltaE window: [-21, 7] MeV applied in ROOT")
alg_pipipi.note(:mbc_signal, "M_BC signal region: [2.282, 2.291] GeV/c^2 applied in ROOT")
alg_pipipi.note(:best_combination, "Best Lambda_c+ candidate chosen by minimum |DeltaE|")
alg_pipipi.note(:peaking_bkg, "No significant peaking background for this mode")
alg_pipipi.note(:deltaE_suppression, "DeltaE(p pi- K+ pi+ pi-) < -40 MeV for bkg suppression")
alg_pipipi.note(:multi_energy, "7 energy points fitted simultaneously with same BF constraint; ARGUS bkg")
alg_pipipi.note(:upper_limit, "Upper limit at 90% CL: B < 6.5 x 10^-4")

alg_pipipi.execute_on(all_data + all_inc_mc + [sig_mc_pipipi])