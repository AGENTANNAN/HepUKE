# Paper: 2109.12660v1
# Title: D_s+ -> pi+ pi0 pi0 amplitude analysis and BF measurement
# Data: 6.32 fb^{-1} at sqrt(s) = 4.178-4.226 GeV (6 energy points)
# Production: e+e- -> D_s*± D_s∓, D_s*+ -> gamma D_s+
#
# Analysis type: Tag-based (tag_side = D_s species)
# Single tag on D_s- via 7 hadronic modes, signal side D_s+ -> pi+ pi0 pi0
# Two-stage: (1) BF measurement via DT/ST, (2) Amplitude analysis (post-DSL)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Energy points (4.178-4.226 GeV) — BOSS 703 datasets
data_4178 = DatasetManager.real_data.find("703_4180")
data_4189 = DatasetManager.real_data.find("703_4190")
data_4199 = DatasetManager.real_data.find("703_4200")
data_4209 = DatasetManager.real_data.find("703_4210")
data_4219 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

all_data = [data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]
all_incMC = [incMC_4178, incMC_4189, incMC_4199, incMC_4209, incMC_4219, incMC_4226]

# ============================================================
# TagAnalysis: D_s+ -> pi+ pi0 pi0
# Tag side: D_s- reconstructed in 7 hadronic modes
# D_s species = :Ds
# Signal side: pi+ pi0 pi0 (1 charged track + 4 photons)
#
# Tag modes (from Table 2 of paper):
#   Mode 1: D_s- -> K_S0 K-
#   Mode 2: D_s- -> K+ K- pi-
#   Mode 3: D_s- -> K_S0 K+ pi0
#   Mode 4: D_s- -> K_S0 K- pi- pi+
#   Mode 5: D_s- -> K_S0 K+ pi- pi-
#   Mode 6: D_s- -> pi- eta' (eta' -> pi+ pi- eta)
#   Mode 7: D_s- -> K- pi+ pi-
#
# Single-tag selection: M_tag in mass window (Table 2), M_rec in range (Table 1)
# M_rec = recoil mass of tag D_s- candidate, used to confirm D_s+ existence
# ============================================================

alg_ds_pipi0pi0 = TagAnalysis.new("DsToPiPi0Pi0")
alg_ds_pipi0pi0.set_header(["DsToPiPi0Pi0Alg/DsToPiPi0Pi0.h"])
                .set_constant({ "ECMS" => [:double, 4.178] })  # representative energy

# Tag side: D_s- reconstructed via 7 hadronic tag modes
alg_ds_pipi0pi0.tag_side(:Ds) do |t|
  t.modes :DstoKSK, :DstoKKPi, :DstoKSKPi0,
          :DstoKSK3Pi, :DstoKSKPiPi,
          :DstoPiEtap, :DstoKPiPi
  t.rank_by :inv   # rank tag candidates by invariant mass
end

# Signal side: D_s+ -> pi+ pi0 pi0
# Production: e+e- -> D_s*± D_s∓ -> gamma D_s+ D_s-
# Transition photon from D_s* decay always present
alg_ds_pipi0pi0.signal_side do |s|
  s.photons 5                        # 4 from 2 pi0 + 1 transition photon (D_s* -> gamma D_s)
  s.charged(pip: 1)                  # pi+ from D_s+ decay
  s.min_photon_energy 0.025
  s.min_photon_angle 10.0
end

# Fit: 4C + pi0 mass constraints + D_s mass constraint
alg_ds_pipi0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_ds_pipi0pi0.note(:tag_selection, "Single tag: M_tag within mode-dependent window (Table 2). " \
                                     "M_rec = recoil mass against tag D_s-; window per Table 1. " \
                                     "Best DT candidate chosen by min chi2 of 8C kinematic fit.")
alg_ds_pipi0pi0.note(:amplitude_analysis, "Amplitude analysis further selects: 8C kinematic fit " \
                                           "(4C + 2 pi0 mass constraints + tag D_s mass + D_s* mass). " \
                                           "Transition photon E_gamma < 0.18 GeV. " \
                                           "M_recoil(gamma_trans, D_s+) in [1.952, 1.995] GeV/c2. " \
                                           "eta veto: reject if any gammagamma pair consistent with eta mass. " \
                                           "D0->Kpipi0 veto: reject if D0 mass from K+pi+pi0 combinations. " \
                                           "K_S0->pi0pi0 veto: M(pi0pi0) not in [0.458, 0.520] GeV/c2. " \
                                           "Signal mass window: [1.925, 1.985] GeV/c2.")
alg_ds_pipi0pi0.note(:bf_method, "BF = N_DT / (N_ST × epsilon_sig). " \
                                  "N_ST from fit to M_tag distributions. " \
                                  "N_DT from counting in M_sig signal region.")
alg_ds_pipi0pi0.note(:transition_photon, "D_s*+ -> gamma D_s+; transition photon always present. " \
                                          "Cross section for D_s*D_s about 20x larger than D_s D_s.")
alg_ds_pipi0pi0.note(:tag_mode_specifics, "Tag mass windows: K_S0K- [1.948,1.991]; K+K-pi- [1.950,1.986]; " \
                                           "K_S0K+pi0 [1.946,1.987]; K_S0K-pi-pi+ [1.958,1.980]; " \
                                           "K_S0K+pi-pi- [1.953,1.983]; pi-eta' [1.940,1.996]; " \
                                           "K-pi+pi- [1.953,1.986] GeV/c2. " \
                                           "Photons: E>25(50) MeV barrel(endcap), |cos_theta|<0.80(0.86-0.92). " \
                                           "pi0 mass: [0.115,0.150] GeV/c2, eta mass: [0.490,0.580] GeV/c2. " \
                                           "Kinematic fit chi2<30 for pi0/eta mass constraints. " \
                                           "eta' mass: [0.946,0.970] GeV/c2.")

alg_ds_pipi0pi0.apply
alg_ds_pipi0pi0.execute_on(all_data + all_incMC)