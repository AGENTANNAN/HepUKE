# BESIII Analysis: e+e- -> K- Xi+ Lambda/Sigma0 cross sections
# Paper: 2406.18183v2 — 35 energy points from sqrt(s)=3.510 to 4.914 GeV
# Method: Partial reconstruction using recoil mass of K- Xi+ system

### Dataset preparation ###

# Real data at scan energy points (BOSS 703, 705, 706, 707, 709, 712)
# Energy points matched to BES3_dataset.md
data_3650 = DatasetManager.real_data.find("709_3650")     # 3.650 GeV
data_3773 = DatasetManager.real_data.find("712_3773")     # 3.773 GeV
data_3872 = DatasetManager.real_data.find("703_3872")     # 3.867 GeV
data_4009 = DatasetManager.real_data.find("703_4009")     # 4.009 GeV
data_4130 = DatasetManager.real_data.find("705_4130")     # 4.130 GeV
data_4160 = DatasetManager.real_data.find("705_4160")     # 4.160 GeV
data_4180 = DatasetManager.real_data.find("703_4180")     # 4.180 GeV
data_4190 = DatasetManager.real_data.find("703_4190")     # 4.190 GeV
data_4200 = DatasetManager.real_data.find("703_4200")     # 4.200 GeV
data_4210 = DatasetManager.real_data.find("703_4210")     # 4.210 GeV
data_4220 = DatasetManager.real_data.find("703_4220")     # 4.220 GeV
data_4230 = DatasetManager.real_data.find("703_4230")     # 4.230 GeV
data_4237 = DatasetManager.real_data.find("703_4237")     # 4.237 GeV
data_4246 = DatasetManager.real_data.find("703_4246")     # 4.246 GeV
data_4260 = DatasetManager.real_data.find("703_4260")     # 4.260 GeV
data_4290 = DatasetManager.real_data.find("705_4290")     # 4.290 GeV
data_4315 = DatasetManager.real_data.find("705_4315")     # 4.315 GeV
data_4340 = DatasetManager.real_data.find("705_4340")     # 4.340 GeV
data_4360 = DatasetManager.real_data.find("703_4360")     # 4.360 GeV
data_4380 = DatasetManager.real_data.find("705_4380")     # 4.380 GeV
data_4400 = DatasetManager.real_data.find("705_4400")     # 4.400 GeV
data_4420 = DatasetManager.real_data.find("703_4420")     # 4.420 GeV
data_4440 = DatasetManager.real_data.find("705_4440")     # 4.440 GeV
data_4600 = DatasetManager.real_data.find("703_4600")     # 4.600 GeV
data_4640 = DatasetManager.real_data.find("706_4640")     # 4.640 GeV
data_4660 = DatasetManager.real_data.find("706_4660")     # 4.660 GeV
data_4680 = DatasetManager.real_data.find("706_4680")     # 4.680 GeV
data_4700 = DatasetManager.real_data.find("706_4700")     # 4.700 GeV
data_4750 = DatasetManager.real_data.find("707_4750")     # 4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")     # 4.780 GeV
data_4914 = DatasetManager.real_data.find("707_4914")     # 4.914 GeV

# Combine all scan points
scan_points = [
  data_3650, data_3773, data_3872, data_4009,
  data_4130, data_4160, data_4180, data_4190,
  data_4200, data_4210, data_4220, data_4230,
  data_4237, data_4246, data_4260, data_4290,
  data_4315, data_4340, data_4360, data_4380,
  data_4400, data_4420, data_4440, data_4600,
  data_4640, data_4660, data_4680, data_4700,
  data_4750, data_4780, data_4914
]

# Decay card for signal process: e+e- -> K- anti-Xi- Lambda
# (Charge-conjugate channel implied; Lambda/Sigma0 inferred from recoil)
# KKMC + psi(4260) as top mother for continuum production
decay_card_KXiLambda = <<~DECAYCARD
  Decay psi(4260)
  1.0 K- anti-Xi- Lambda PHSP;
  Enddecay

  Decay anti-Xi-
  1.0 anti-Lambda pi+ PHSP;
  Enddecay

  Decay anti-Lambda
  1.0 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal at all scan points
exMC_KXiLambda = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_KXiLambda"
  config.events        = 200_000
  config.decay_card    = decay_card_KXiLambda
  config.cross_section = :default
end

# Its inclusive MC counterparts (one per energy point; listing representative)
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
# ... (inclusive MC for all scan points)

### Event selection (BOSS) ###

alg_KXiLambda = Algorithm.new("KXiLambdaRecoil")
alg_KXiLambda.set_header(["KXiLambdaRecoilAlg/KXiLambdaRecoil.h"])

# The analysis is a multi-energy scan; ECMS is set per energy point
# by execute_on automatically (no single set_constant "ECMS" needed)

event_selection = Selection.new

# Charged track selection: |cos(theta)| < 0.93, >=2 positive, >=2 negative tracks
event_selection.select_track do
  cos_theta   0.93
  nChrp       ">=2"
  nChrn       ">=2"
end

# Particle identification: likelihood-based dE/dx + TOF
# Proton: L(p) > L(K) and L(p) > L(pi)
# Kaon:   L(K) > L(pi)
# Pion:   L(pi) > L(K)
# Require >=1 anti-proton, >=1 K-, >=2 pi+
event_selection.pid(method: :probability) do
  prob_cut   0.001
  identify :proton, against: [:kaon, :pion]
  identify :kaon,   against: [:pion]
  identify :pion,   against: [:kaon]
  nprm       ">=1"    # at least one anti-proton (from anti-Lambda decay)
  nkm        ">=1"    # at least one K-
  npip       ">=2"   # at least two pi+ (one from anti-Lambda, one from anti-Xi-)
end

# Reconstruct anti-Lambda -> anti-p pi+ via secondary vertex fit
# Best candidate selected by minimizing mass difference wrt PDG Lambda mass
event_selection.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Partial reconstruction: reconstruct K- + anti-Xi- system,
# infer Lambda (recID 3) from recoil mass
# RecID mapping from decay card:
#   0: psi(4260) [top mother, replaced by P4_cms]
#   1: K-
#   2: anti-Xi-  [intermediate; daughters expanded automatically]
#   3: Lambda    [MISSED — inferred from recoil]
#   4: anti-Lambda [pre-built by secondary_vertex_fit as Lambda_bar]
#   5: pi+       [from anti-Xi- decay]
#   6: anti-p-   [from anti-Lambda decay]
#   7: pi+       [from anti-Lambda decay]
event_selection.partial_miss([3]) do
  # Recoil mass window: 1.0 – 1.3 GeV/c² covers both Lambda and Sigma0 peaks
  require_recoil_mass 1.0, 1.3
end

alg_KXiLambda.with_decay_card(decay_card_KXiLambda).apply(event_selection)

# Notes for inexpressible BOSS-side procedures
alg_KXiLambda
  .note(:decay_length_cut,
    "Decay length of reconstructed Lambda_bar and anti-Xi- candidates required > 0 " \
    "to suppress non-Lambda and non-Xi background; applied via secondary vertex fit " \
    "decay length significance in BOSS")
  .note(:lambda_mass_window,
    "anti-p pi+ invariant mass window: +/- 8 MeV/c2 around nominal Lambda mass, " \
    "optimized by FOM S/sqrt(S+B); applied in BOSS after secondary vertex fit")
  .note(:xi_mass_window,
    "pi+ anti-Lambda invariant mass window: +/- 6 MeV/c2 around nominal anti-Xi- mass, " \
    "optimized by FOM; applied in BOSS after Xi reconstruction")
  .note(:combined_mass_minimization,
    "Best Xi and Lambda candidates chosen by minimizing combined mass difference " \
    "|M(pi+ anti-Lambda) - m(anti-Xi-)| + |M(pi+ anti-p) - m(anti-Lambda)|; " \
    "DSL uses independent by_minimizing_mass_difference per vertex fit")
  .note(:conjugate_mode,
    "Charge-conjugate channel e+e- -> K+ Xi- anti-Lambda is always implied; " \
    "factor 2 in cross section formula accounts for both charge modes")
  .note(:sigma0_channel,
    "e+e- -> K- anti-Xi- Sigma0 also measured; Sigma0 -> Lambda gamma. " \
    "Sigma0 peak appears in the same recoil mass spectrum at higher mass. " \
    "Same selection serves both Lambda and Sigma0 channels")
  .note(:isr_correction,
    "ISR correction factor (1+delta) obtained via QED calculation [Kuraev-Fadin]; " \
    "iterative procedure reweights signal MC until (1+delta)*epsilon converges <0.1%")

# Execute algorithm on the scan datasets + inclusive MC + exclusive MC
all_datasets = scan_points + [exMC_KXiLambda].flatten
alg_KXiLambda.execute_on(all_datasets)