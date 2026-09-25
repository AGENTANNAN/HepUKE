# Paper: 2309.05484v2 — First observation of Lambda_c+ -> Sigma- K+ pi+
# Method: Double-tag (DT) with 10 ST hadronic modes for Lambda_c-
# Data: 4.5 fb-1 at 7 energy points 4.600-4.699 GeV
# Result: Observation at 5.4 sigma; BFs measured for Sigma- K+ pi+, Xi- K+ pi+, and Xi*0 K+

### Dataset preparation ###
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 7 energy points spanning 4.600-4.699 GeV
data_4600 = DatasetManager.real_data.find("706_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_4600 = DatasetManager.inclusive_mc.find("706_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

all_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card: e+e- -> Lambda_c+ Lambda_c- via KKMC
# Signal: Lambda_c+ -> Sigma- K+ pi+, Sigma- -> n pi- (neutron inferred from missing energy)
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 anti-Lambda_c- Lambda_c+ PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 Sigma- K+ pi+ PHSP;
    Enddecay

    Decay Sigma-
    1.000 n0 pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "Lambdac_Sigma_K_pi"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# 10 tag modes for Lambda_c- (from paper Table I)
# Mapped to authoritative DTagAlg channel names:
#   pK+pi-     -> LambdacPtoKPiP
#   pKS0       -> LambdacPtoKsP
#   pK+pi-pi0  -> LambdacPtoKPiPi0P
#   pKS0pi0    -> LambdacPtoKsPi0P
#   pKS0pi+pi- -> LambdacPtoKsPiPiP
#   Lambda pi- -> LambdacPtoLambdaPi
#   Lambda pi-pi0 -> LambdacPtoLambdaPiPi0
#   Lambda pi+pi-pi- -> LambdacPtoLambdaPiPiPi
#   Sigma0 pi- -> LambdacPtoPiSIGMA0LambdaGam
#   Sigma- pi+pi- -> NOT IN VOCABULARY (dropped)
lc_tag_modes = [:LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoKPiPi0P,
                :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP, :LambdacPtoLambdaPi,
                :LambdacPtoLambdaPiPi0, :LambdacPtoLambdaPiPiPi,
                :LambdacPtoPiSIGMA0LambdaGam]

### Algorithm: Lambda_c+ -> Sigma- K+ pi+ (Sigma- -> n pi-) ###
alg = TagAnalysis.new("LambdacSigmaKPi")
alg.set_header(["LambdacSigmaKPiAlg/LambdacSigmaKPi.h"])
alg.set_constant({ "ECMS" => [:double, 4.650] })

# Tag side: Lambda_c- reconstructed via 10 hadronic modes
alg.tag_side(:"Lambda_c+") do |t|
  t.modes(*lc_tag_modes)
  t.charm -1       # tag the Lambda_c- (charm = -1)
end

# Signal side: Lambda_c+ -> Sigma- K+ pi+
# Sigma- decays via Sigma- -> n pi- (neutron missing)
# Charged: K+, pi+ (from Lambda_c+ decay), pi- (from Sigma- decay)
alg.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1   # K+(+1) + pi+(+1) + pi-(-1) = +1
  s.missing :n0        # neutron from Sigma- -> n pi- (massive missing)
end

# Kinematic fit: 4C energy-momentum conservation + Lambda_c mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg
  .note(:tag_mode_unavailable,
    "Sigma- pi+pi- tag mode (Lambda_c- -> Sigma- pi+pi-) not found in BOSS 7.0.6/7.1.2 DTagAlg vocabulary; dropped from ST selection")
  .note(:st_yields,
    "ST Lambda_c- yields from previous BESIII measurements: total 105249 +/- 386 events across all 10 tag modes and 7 energy points")
  .note(:track_selection,
    "tracks required within |cos_theta| < 0.93; K+/pi+ tracks: |Vz| < 10 cm, |Vr| < 1 cm; pi- from Sigma-: |Vz| < 20 cm")
  .note(:pid,
    "PID uses combined dE/dx and TOF; K+ and pi+ identified by PID probability")
  .note(:vertex_fit,
    "vertex fit performed to K+ and pi+ candidates; momenta updated by fit used in subsequent analysis")
  .note(:extra_track_veto,
    "no extra charged tracks with |cos_theta| < 0.93, |Vr| < 1 cm, |Vz| < 20 cm")
  .note(:mrec_lambdac_window,
    "M_rec(Lambda_c-) required within (2.275, 2.310) GeV/c^2 to suppress qqbar continuum background")
  .note(:sigma_plus_background_veto,
    "events with M_rec(H+) in (1.15, 1.24) GeV/c^2 vetoed to suppress Lambda_c+ -> Sigma+ K+ pi- background")
  .note(:ks_veto,
    "events with M(pi+ pi-) in (0.48, 0.52) GeV/c^2 vetoed to suppress Lambda_c+ -> n KS0 K+ background")
  .note(:mrec_sigma_cut,
    "M_rec(H-) > 1.15 GeV/c^2 to suppress qqbar and non-signal Lambda_c+ Lambda_c- backgrounds")
  .note(:signal_extraction,
    "signal yields from 2D unbinned maximum-likelihood fit to M_rec(H-) vs M_rec(B0); " \
    "Sigma- and neutron signals identified in 2D distribution; " \
    "also measures Lambda_c+ -> Xi- K+ pi+ (128 +/- 13 events) and Lambda_c+ -> Xi*0 K+ (54 +/- 8 events)")
  .note(:bf_calculation,
    "BF = N_obs / sum(N_ST_ij * epsilon_DT_ij / epsilon_ST_ij); " \
    "ST efficiencies cancel in DT method; " \
    "ST yields and efficiencies from Refs.[31,32]")
  .note(:mc_reweighting,
    "key distributions of ST modes reweighted to match data; " \
    "Lambda_c+ -> Xi- K+ pi+ signal MC reweighted to match data + Lambda_c+ -> Xi(1530)0 K+ component")
  .note(:inclusive_mc_scale,
    "inclusive MC normalized with scale factor 0.034 from M_BC sideband comparison (2.10, 2.25) GeV/c^2")

alg.with_decay_card(decay_card).apply
alg.execute_on(all_data + all_incMC + exMC_signal)