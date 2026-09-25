# Paper: 2212.03753v2 — Improved measurement of absolute BF of inclusive semileptonic Lambda_c+ decay
# Process: Lambda_c+ → X e+ nu_e (inclusive semileptonic)
# TagAnalysis: ST Lambda_c- (hadronic tags) + signal side e+ + missing nu_e
# 7 energy points: 4.600, 4.612, 4.628, 4.640, 4.661, 4.682, 4.698 GeV

### Dataset preparation ###
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4610")  # ~4.612
data_4628 = DatasetManager.real_data.find("706_4620")  # ~4.628
data_4640 = DatasetManager.real_data.find("706_4640")  # ~4.641
data_4660 = DatasetManager.real_data.find("706_4660")  # ~4.661
data_4680 = DatasetManager.real_data.find("706_4680")  # ~4.682
data_4700 = DatasetManager.real_data.find("706_4700")  # ~4.699

scan_data = [data_4600, data_4612, data_4628, data_4640, data_4660, data_4680, data_4700]

scan_incMC = scan_data.map do |ds|
  DatasetManager.inclusive_mc.find(ds.sample_name)
end

# Decay card for signal MC: Lambda_c+ → X e+ nu_e (inclusive)
decay_card = <<~DECAYCARD
    Decay Lambda_c+
    1.0 X e+ nu_e PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_Lambdac_Xenu"
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ============================================================
# TagAnalysis: ST Lambda_c- → hadronic modes, signal side e+ + missing nu_e
# ST + missing pattern (semileptonic)
# ============================================================
alg_lambdac_sl = TagAnalysis.new("LambdacSL")
alg_lambdac_sl.set_header(["LambdacSLAlg/LambdacSL.h"])
                .set_constant({ "ECMS" => [:double, 4.640] })
                .note(:tag_mode_unavailable,
                  "Paper uses 12 specific Lambda_c- hadronic tag modes:
                   pbar K_S^0, pbar K+ pi-, pbar K_S^0 pi0, pbar K+ pi- pi0,
                   pbar K_S^0 pi+ pi-, anti-Lambda pi-, anti-Lambda pi- pi0,
                   anti-Lambda pi- pi+ pi-, anti-Sigma0 pi-, anti-Sigma- pi0,
                   anti-Sigma- pi+ pi-, Sigma0 pi- pi0.
                   Using mode_group :hadronic as closest approximation;
                   unavailable modes noted.")
                .note(:signal_side_description,
                  "Signal side: one positron (e+) identified via E/p > 0.8 and
                   PID likelihood cuts (L'_e > 0.001, L'_e/(L'_e+L'_π+L'_K) > 0.8).
                   Missing neutrino nu_e (massless). Positron momentum bins
                   50 MeV/c from 200-1000 MeV/c. WS subtraction for charge-symmetric
                   backgrounds. PID unfolding and tracking efficiency unfolding
                   applied in ROOT analysis stage, not in BOSS selection.")

alg_lambdac_sl.tag_side(:Lambdac) do |t|
  # 12 hadronic tag modes for Lambda_c- (anti-charm baryon)
  # Include all available hadronic Lambda_c modes
  t.mode_group :hadronic
  # Explicitly include known modes; mode_group :hadronic covers the hadronic block
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP
  # Pin to anti-charm (Lambda_c-)
  t.charm(-1)
end

alg_lambdac_sl.signal_side do |s|
  # One positron on signal side
  s.charged(ep: 1, at_least: true)
  s.require_charge 1
  # Missing massless neutrino
  s.missing :nu_e
end

alg_lambdac_sl.fit do |f|
  # 4C: tag(Lambda_c-) + e+ + nu_e = ecms_lab
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_lambdac_sl.apply

all_datasets = scan_data + scan_incMC + exMC_signal
alg_lambdac_sl.execute_on(all_datasets)