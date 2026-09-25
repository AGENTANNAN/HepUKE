# ============================================================
# Paper: Search for e+e- → ηc η π+π- at BESIII
# arXiv: 2011.13850v2
# ============================================================
# ηc reconstructed in 16 exclusive decay modes.
# 16 independent Algorithm objects, one per ηc decay mode.
# Representative modes are given here; the pattern is identical
# for the remaining modes with varying track/photon counts.

###
### Datasets: 5 energy points, BOSS 703
###
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

data_points = [data_4230, data_4260, data_4360, data_4420, data_4600]
inc_mc_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# ============================================================
# MODE 1: ηc → π+ π- K+ K-
# Final state: π+π- (from e+e-) + γγ (η→γγ) + π+π- + K+K- (ηc decay)
# Tracks: 2 pos + 2 neg from η+ρ → 4 pions; 4 from ηc → 8 total = 4 pos + 4 neg
# Photons: 2 (from η→γγ)
# ============================================================
decay_card_m1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- K+ K-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_pipi_KK_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m1
  config.cross_section = :default
end

alg_m1 = Algorithm.new("EtacEtaPipPim_PiPiKK")
alg_m1.set_header(["EtacEtaPipPim_PiPiKKAlg/EtacEtaPipPim_PiPiKK.h"])
alg_m1.set_constant({})

sel_m1 = Selection.new
sel_m1.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       "==4"
  nChrn       "==4"
  nNet        "==0"
end
sel_m1.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              "==2"
end
# PID for kaons
sel_m1.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==2"
  nkm "==2"
end
sel_m1.remove([:kp <= :chrgp, :km <= :chrgn])
# Remaining tracks are pions
sel_m1.assign({chrgp: :pip, chrgn: :pim})
# 6C fit: 4C + η mass constraint
sel_m1.kinematic_fit([:pip, :pim, :pip, :pim, :kp, :km, :gamma, :gamma]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
end
alg_m1.note(:ks_vertex, "K_S0 vertex fit L/sigma_L > 2; mass window ±15 MeV/c2 around K_S0 nominal")
alg_m1.note(:pi0_eta_reconstruction, "π0 mass window [110,150] MeV/c2; η mass window [500,570] MeV/c2; mass-constrained via Kalman fit")
alg_m1.note(:chi2_selection, "Kinematic fit χ2 cut chosen per final state to retain 90% of signal events")
alg_m1.note(:multiple_candidates, "Multiple candidates per event with same χ2 kept; ROOT-level fit resolves via smooth background")
alg_m1.with_decay_card(decay_card_m1).apply(sel_m1)
alg_m1.execute_on(data_points + inc_mc_points + exMC_m1)

# ============================================================
# MODE 2: ηc → 2(K+ K-)
# Final state: π+π- + γγ + 2(K+K-) → 4 pos + 4 neg + 2 photons
# ============================================================
decay_card_m2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 K+ K- K+ K-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_4K_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m2
  config.cross_section = :default
end

alg_m2 = Algorithm.new("EtacEtaPipPim_4K")
alg_m2.set_header(["EtacEtaPipPim_4KAlg/EtacEtaPipPim_4K.h"])
alg_m2.set_constant({})
sel_m2 = Selection.new
sel_m2.select_track do
  cos_theta   0.93; Vz 10.0; Vr 1.0; nChrp "==4"; nChrn "==4"; nNet "==0"
end
sel_m2.select_photon do
  tdc_emc_start 0; tdc_emc_end 14; energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0; nGam "==2"
end
sel_m2.pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; nkp "==4"; nkm "==4" }
sel_m2.remove([:kp <= :chrgp, :km <= :chrgn])
sel_m2.assign({chrgp: :pip, chrgn: :pim})
sel_m2.kinematic_fit([:pip, :pim, :pip, :pim, :kp, :km, :kp, :km, :gamma, :gamma]) do
  nominal; constrain_four_momentum; invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta); chi2_cut 200
end
alg_m2.with_decay_card(decay_card_m2).apply(sel_m2)
alg_m2.execute_on(data_points + inc_mc_points + exMC_m2)

# ============================================================
# MODE 3: ηc → 2(π+ π-)
# Final state: π+π- + γγ + 2(π+π-) → 4 pos + 4 neg + 2 photons
# ============================================================
decay_card_m3 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- pi+ pi-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m3 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_4pi_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m3
  config.cross_section = :default
end

alg_m3 = Algorithm.new("EtacEtaPipPim_4Pi")
alg_m3.set_header(["EtacEtaPipPim_4PiAlg/EtacEtaPipPim_4Pi.h"])
alg_m3.set_constant({})
sel_m3 = Selection.new
sel_m3.select_track do
  cos_theta   0.93; Vz 10.0; Vr 1.0; nChrp "==4"; nChrn "==4"; nNet "==0"
end
sel_m3.select_photon do
  tdc_emc_start 0; tdc_emc_end 14; energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0; nGam "==2"
end
sel_m3.assign({chrgp: :pip, chrgn: :pim})
sel_m3.kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim, :gamma, :gamma]) do
  nominal; constrain_four_momentum; invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta); chi2_cut 200
end
alg_m3.with_decay_card(decay_card_m3).apply(sel_m3)
alg_m3.execute_on(data_points + inc_mc_points + exMC_m3)

# ============================================================
# MODE 4: ηc → K+ K- π0
# Final state: π+π- + γγ(η) + K+K- + γγ(π0)
# Tracks: 2 pos + 2 neg (pions from e+e-) + 2 kaons → 3 pos + 3 neg
# Photons: 4 (2 from η, 2 from π0)
# ============================================================
decay_card_m4 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 K+ K- pi0  PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m4 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_KKpi0_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m4
  config.cross_section = :default
end

alg_m4 = Algorithm.new("EtacEtaPipPim_KKPi0")
alg_m4.set_header(["EtacEtaPipPim_KKPi0Alg/EtacEtaPipPim_KKPi0.h"])
alg_m4.set_constant({})
sel_m4 = Selection.new
sel_m4.select_track do
  cos_theta   0.93; Vz 10.0; Vr 1.0; nChrp "==3"; nChrn "==3"; nNet "==0"
end
sel_m4.select_photon do
  tdc_emc_start 0; tdc_emc_end 14; energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0; nGam ">=4"
end
sel_m4.pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; nkp "==1"; nkm "==1" }
sel_m4.remove([:kp <= :chrgp, :km <= :chrgn])
sel_m4.assign({chrgp: :pip, chrgn: :pim})
sel_m4.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 30
  npi0 ">=1"
end
sel_m4.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 30
  neta ">=1"
end
sel_m4.kinematic_fit([:pip, :pim, :pip, :pim, :kp, :km, :eta, :pi0]) do
  nominal; constrain_four_momentum; chi2_cut 200
end
alg_m4.with_decay_card(decay_card_m4).apply(sel_m4)
alg_m4.execute_on(data_points + inc_mc_points + exMC_m4)

# ============================================================
# MODE 5: ηc → K+ K- η
# Final state: π+π- + γγ(η) + K+K- + γγ(η_c daughter η)
# Tracks: 2 pos + 2 neg (pions) + 2 kaons → 3 pos + 3 neg
# Photons: 4 (2 from primary η, 2 from secondary η)
# ============================================================
decay_card_m5 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 K+ K- eta  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m5 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_KKeta_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m5
  config.cross_section = :default
end

alg_m5 = Algorithm.new("EtacEtaPipPim_KKeta")
alg_m5.set_header(["EtacEtaPipPim_KKetaAlg/EtacEtaPipPim_KKeta.h"])
alg_m5.set_constant({})
sel_m5 = Selection.new
sel_m5.select_track do
  cos_theta   0.93; Vz 10.0; Vr 1.0; nChrp "==3"; nChrn "==3"; nNet "==0"
end
sel_m5.select_photon do
  tdc_emc_start 0; tdc_emc_end 14; energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0; nGam ">=4"
end
sel_m5.pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; nkp "==1"; nkm "==1" }
sel_m5.remove([:kp <= :chrgp, :km <= :chrgn])
sel_m5.assign({chrgp: :pip, chrgn: :pim})
sel_m5.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 30
  neta ">=2"
end
sel_m5.kinematic_fit([:pip, :pim, :pip, :pim, :kp, :km, :eta, :eta]) do
  nominal; constrain_four_momentum; chi2_cut 200
end
alg_m5.with_decay_card(decay_card_m5).apply(sel_m5)
alg_m5.execute_on(data_points + inc_mc_points + exMC_m5)

# ============================================================
# MODE 6: ηc → p pbar
# Final state: π+π- + γγ(η) + p + pbar
# Tracks: 2 pions + p + pbar → 2 pos (π+, p) + 2 neg (π-, pbar)
# Photons: 2 (from η→γγ)
# ============================================================
decay_card_m6 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta_c eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    Decay eta_c
    1.000 p+ anti-p-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_m6 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etac_ppbar_etapp"
  config.events        = 200_000
  config.decay_card    = decay_card_m6
  config.cross_section = :default
end

alg_m6 = Algorithm.new("EtacEtaPipPim_ppbar")
alg_m6.set_header(["EtacEtaPipPim_ppbarAlg/EtacEtaPipPim_ppbar.h"])
alg_m6.set_constant({})
sel_m6 = Selection.new
sel_m6.select_track do
  cos_theta   0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0"
end
sel_m6.select_photon do
  tdc_emc_start 0; tdc_emc_end 14; energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0; nGam "==2"
end
# Identify protons; remaining are pions
sel_m6.pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
sel_m6.remove([:prp <= :chrgp, :prm <= :chrgn])
sel_m6.assign({chrgp: :pip, chrgn: :pim})
sel_m6.kinematic_fit([:pip, :pim, :prp, :prm, :gamma, :gamma]) do
  nominal; constrain_four_momentum; invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta); chi2_cut 200
end
alg_m6.with_decay_card(decay_card_m6).apply(sel_m6)
alg_m6.execute_on(data_points + inc_mc_points + exMC_m6)

# ============================================================
# Remaining 10 modes follow identical patterns; omitted for brevity.
# Modes: 3(π+π-), K_S0 K± π∓, K_S0 K± π∓ π+π-, π+π-η,
#        π+π-π0π0, 2(π+π-)η, 2(π+π-π0), K+K-2(π+π-),
#        pp̄π0, pp̄π+π-
# ============================================================