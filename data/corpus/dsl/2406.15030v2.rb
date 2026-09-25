# Paper: 2406.15030v2
# Search for e+e- → φχc1(3872) at 4.914 and 4.946 GeV
# φ → K+K-, χc1(3872) → ρ0J/ψ → π+π-ℓ+ℓ-, ℓ=e/μ
# Two event classes: 6-track (no missing) and 5-track (missing K)

# Datasets: BOSS 707, two energy points
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

data_points = [data_4914, data_4946]

# Decay card for e+e- → φ χc1(3872), φ → K+K-, χc1(3872) → ρ0 J/ψ → π+π- ℓ+ℓ-
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 phi chi_c1_3872 PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay chi_c1_3872
    1.0000 rho0 J/psi PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- VSS;
    Enddecay

    Decay J/psi
    1.0000 l+ l- PHSP;
    Enddecay

    End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "ee_phi_chi_c1_3872"
  config.events = 50000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ===============================
# Algorithm 1: 6-track events
# ===============================
alg_6trk = Algorithm.new("PhiChic13872_6Trk")
alg_6trk.set_header(["PhiChic13872_6TrkAlg/PhiChic13872_6Trk.h"])
        .set_constant({ "ECMS" => [:double, 4.918] })
        .note(:ecms_per_dataset,
          "ECMS is 4918.02 MeV for 707_4914 and 4950.93 MeV for 707_4946; per-run MeasuredEcmsSvc handles the actual value")

# 6-track selection: K+ K- π+ π- ℓ+ ℓ-, zero net charge
sel_6trk = Selection.new
sel_6trk.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==3"
  nChrn "==3"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
end
# Lepton identification: high-momentum tracks treated as leptons
# Electron/muon separation via EMC energy (E>0.8 GeV → e, E<0.4 GeV → μ)
# Note: the paper uses custom thresholds (0.8/0.4) vs DSL default (0.6)
.identify_high_momentum_leptons do
  treat_as_lepton_if_momentum_above 1.0
  treat_as_electron_if_energy_above 0.6
end
# 4C kinematic fit: K+ K- π+ π- ℓ+ ℓ- under 4-momentum conservation
.kinematic_fit([:kp, :km, :pip, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 150
end

alg_6trk
  .note(:lepton_id_custom,
    "Paper uses EMC E>0.8 GeV for electron, E<0.4 GeV for muon, and MUC depth>30cm for muon.
     DSL identify_high_momentum_leptons uses a single E/eraw>0.6 threshold.
     The custom PID is implemented in ROOT analysis stage.")
  .note(:mass_windows,
    "J/ψ mass window |M(ℓ+ℓ-) - m(J/ψ)| < 0.04 GeV/c²;
     φ mass window |M(K+K-) - m(φ)| < 0.015 GeV/c²;
     both applied in ROOT after kinematic fit")
  .with_decay_card(decay_card)
  .apply(sel_6trk)

alg_6trk.execute_on(data_points + [incMC_4914, incMC_4946] + exMCs)

# ===============================
# Algorithm 2: 5-track events
# ===============================
alg_5trk = Algorithm.new("PhiChic13872_5Trk")
alg_5trk.set_header(["PhiChic13872_5TrkAlg/PhiChic13872_5Trk.h"])
        .set_constant({ "ECMS" => [:double, 4.918] })
        .note(:ecms_per_dataset,
          "ECMS is 4918.02 MeV for 707_4914 and 4950.93 MeV for 707_4946; per-run MeasuredEcmsSvc handles the actual value")

# 5-track selection: one kaon missing
# Accept K±π∓π+π-ℓ+ℓ- (5 charged tracks, zero net charge)
sel_5trk = Selection.new
sel_5trk.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
  nTot "==5"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
end
.identify_high_momentum_leptons do
  treat_as_lepton_if_momentum_above 1.0
  treat_as_electron_if_energy_above 0.6
end
# 1C kinematic fit with missing kaon
.kinematic_fit([:kp, :pip, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  miss_track_of :kp
  chi2_cut 20
end

alg_5trk
  .note(:lepton_id_custom,
    "Paper uses EMC E>0.8 GeV for electron, E<0.4 GeV for muon, and MUC depth>30cm for muon.
     DSL identify_high_momentum_leptons uses a single E/eraw>0.6 threshold.
     The custom PID is implemented in ROOT analysis stage.")
  .note(:mass_windows,
    "J/ψ mass window |M(ℓ+ℓ-) - m(J/ψ)| < 0.04 GeV/c²;
     φ mass window |M(K+K-) - m(φ)| < 0.015 GeV/c²;
     both applied in ROOT after kinematic fit")
  .note(:missing_k,
    "5-track event class: one kaon missing. Both K+ and K- hypotheses tested;
     candidate with best χ²_1C selected in ROOT stage")
  .with_decay_card(decay_card)
  .apply(sel_5trk)

alg_5trk.execute_on(data_points + [incMC_4914, incMC_4946] + exMCs)