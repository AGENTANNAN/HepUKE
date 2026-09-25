# Cross section measurement of e+e- → f1(1285)π+π- at 45 CM energies from 3.808-4.951 GeV
# f1(1285) → π+π-η, η → γγ
# arXiv:2501.14206v1
# ConExc generator for signal MC (ISR continuum process)

# Representative energy points from the 45-point scan (BOSS versions 703, 705, 706, 707, 712)
# Full list includes 45 points; shown are representative samples spanning the energy range
data_3810 = DatasetManager.real_data.find("703_3810")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4946 = DatasetManager.real_data.find("707_4946")

incMC_3810 = DatasetManager.inclusive_mc.find("703_3810")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

# ConExc card: mode 41 = f1(1285)π+π-
# Omit Particle vpho line — DSL auto-injects per-point CMS energy
decay_card = <<~DECAYCARD
  Decay vpho
  1 ConExc 41;
  Enddecay

  Decay vhdr
  1 f_1 pi+ pi- PHSP;
  Enddecay

  Decay f_1
  1 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Multi-energy exclusive MC via create_exclusive_mc_for
scan_data_points = [data_3810, data_4180, data_4260, data_4420, data_4600, data_4680, data_4946]

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name = "conexc_f1pipi"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm: e+e- → f1(1285)π+π-, f1(1285) → π+π-η, η → γγ
# Final state: 2(π+π-) + γγ = 4 charged pions + 2 photons
alg = Algorithm.new("F1PiPiCrossSection")
alg.set_header(["F1PiPiCrossSectionAlg/F1PiPiCrossSection.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })  # representative; varies per scan point

# ── Selection ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==4"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  # PID: all four charged tracks identified as pions
  .pid(method: :probability) {
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  # Reconstruct η → γγ via Kalman kinematic fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # 4C kinematic fit: π+π-π+π- η
  .kinematic_fit([:pip, :pim, :pip, :pim, :eta]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg.with_decay_card(decay_card).apply(event_selection)

# Additional cuts applied in ROOT:
# - χ²_4C < 100 (tight χ² cut applied in ROOT, loose 200 in BOSS per Rule T3)
# - η mass window: |M(γγ) - M(η)| < 25 MeV/c²
# - If >2 photon candidates, select combination with smallest χ²_4C
# - Signal extraction: fit to M(π+π-η) distribution at each energy point
# - ISR correction factor (1+δ_γ) and VP correction (1+δ_v) from iterative method
# - Born/dressed cross sections computed from σ = N_sig / (L·ε·(1+δ_γ)·(1+δ_v)·B)
alg.note(:root_cuts, "ROOT-level: χ²_4C<100, |M(γγ)-M(η)|<25 MeV/c², best photon combination by min χ²_4C, M(π+π-η) signal fit at each √s")
alg.note(:cross_section, "Born/dressed cross sections via iterative method; ISR factor 1+δ_γ and VP factor 1+δ_v from kkmc/ConExc; 45 energy points 3.808-4.951 GeV")
alg.note(:conexc_usage, "Signal MC generated with ConExc generator (mode 41) for ISR continuum e+e-→f1(1285)π+π-; cross-checked against KKMC, difference <0.2%")

alg.execute_on([data_3810, data_4180, data_4260, data_4420, data_4600, data_4680, data_4946,
                incMC_3810, incMC_4180, incMC_4260, incMC_4420, incMC_4600, incMC_4680, incMC_4946,
                exMC_signal])