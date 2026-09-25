# Paper: arXiv:2401.09225v1
# "First measurements of the absolute branching fraction of Lambda_c(2625)+ -> Lambda_c+ pi+ pi-"
# BESIII: e+e- collision data at sqrt(s) = 4.918 and 4.950 GeV (368.48 pb-1 total)
# Method: TagAnalysis with partial reconstruction (missing anti-Lambda_c-)
# Tag modes: pK-pi+ (LambdacPtoKPiP), pK_S0 (LambdacPtoKsP), Lambda pi+ (LambdacPtoLambdaPi)
# Signal: Lambda_c(2625)+ -> Lambda_c+ pi+ pi- with pi+ pi- reconstructed; anti-Lambda_c- missing

### Dataset loading ###
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4918 = DatasetManager.real_data.find("707_4914")   # sqrt(s) = 4.918 GeV, 219.78 pb-1
data_4950 = DatasetManager.real_data.find("707_4946")   # sqrt(s) = 4.950 GeV, 148.70 pb-1
all_data = [data_4918, data_4950]

incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4950 = DatasetManager.inclusive_mc.find("707_4946")
all_incMC = [incMC_4918, incMC_4950]

### Decay card for signal MC ###
# Production: e+e- -> anti-Lambda_c- Lambda_c(2625)+ at sqrt(s) = 4.918/4.950 GeV
# Signal decay: Lambda_c(2625)+ -> Lambda_c+ pi+ pi- (PHSP)
# Lambda_c+ -> p+ K- pi+ (representative tag mode)
# anti-Lambda_c- -> anti-p- K+ pi- (arbitrary decay; cancels in efficiency)
decay_card = <<~DECAYCARD
  Decay vpho
  1.000 anti-Lambda_c- Lambda_c(2625)+ PHSP;
  Enddecay
  Decay Lambda_c(2625)+
  1.000 Lambda_c+ pi+ pi- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 p+ K- pi+ PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Charge conjugate: e+e- -> Lambda_c+ anti-Lambda_c(2625)-
# anti-Lambda_c(2625)- -> anti-Lambda_c- pi+ pi- (PHSP)
decay_card_cc = <<~DECAYCARD
  Decay vpho
  1.000 Lambda_c+ anti-Lambda_c(2625)- PHSP;
  Enddecay
  Decay anti-Lambda_c(2625)-
  1.000 anti-Lambda_c- pi+ pi- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 p+ K- pi+ PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

### Exclusive MC samples ###
exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc2625_Lc_pipi"
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

exMC_cc = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc2625_cc_Lc_pipi"
  config.events = 500_000
  config.decay_card = decay_card_cc
  config.cross_section = :default
end

### TagAnalysis: Lambda_c(2625)+ -> Lambda_c+ pi+ pi- ###
alg = TagAnalysis.new("Lambdac2625")
alg.set_header(["Lambdac2625Alg/Lambdac2625.h"])
   .set_constant({ "ECMS" => [:double, 4.918] })

# Tag side: reconstruct Lambda_c+ via 3 hadronic single-tag modes
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPi
  t.charm 1
end

# Signal side: pi+ pi- from Lambda_c(2625)+ -> Lambda_c+ pi+ pi- decay
# anti-Lambda_c- is unreconstructed (missing)
alg.signal_side do |s|
  s.photons 0
  s.charged(pip: 1, pim: 1)
  s.require_charge 0
  s.missing :X0, mass: 2.28646    # anti-Lambda_c- (PDG nominal mass, unreconstructed)
end

# 4C kinematic fit: constrain total 4-momentum to ECMS
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

### Analysis notes ###
alg.note(:two_energy_points,
  "Analysis uses e+e- collision data at two c.m. energies: " \
  "sqrt(s) = 4.918 GeV (dataset 707_4914, 219.78 pb-1) and " \
  "4.950 GeV (dataset 707_4946, 148.70 pb-1), for a total integrated " \
  "luminosity of 368.48 pb-1. ECMS constant set to 4.918 GeV as a " \
  "placeholder; the actual per-run beam energy is handled by " \
  "MeasuredEcmsSvc in the BOSS framework. The two energy points are " \
  "analysed separately and the results combined via weighted average.")

alg.note(:partial_reconstruction,
  "Partial reconstruction method: the anti-Lambda_c- baryon is not " \
  "reconstructed. Its 4-momentum is inferred from the known initial " \
  "state and the reconstructed Lambda_c(2625)+ candidate via missing " \
  "mass. The Lambda_c+ is fully reconstructed on the tag side via " \
  "hadronic tag modes. The signal-side pi+ pi- pair (from Lambda_c(2625)+ " \
  "decay) is combined with the tagged Lambda_c+ to form the " \
  "Lambda_c(2625)+ candidate. The missing mass of the anti-Lambda_c- " \
  "is set to 2.28646 GeV/c^2 (PDG nominal mass).")

alg.note(:recoil_mass_variable,
  "M_recoil: the recoil mass of the system recoiling against the " \
  "reconstructed Lambda_c+ pi+ pi- combination, defined as: " \
  "M_recoil = sqrt[ (E_cms - E_Lc+ - E_pi+_sig - E_pi-_sig)^2 - " \
  "|p_cms - p_Lc+ - p_pi+_sig - p_pi-_sig|^2 ]. " \
  "This variable peaks at the anti-Lambda_c- mass for signal events. " \
  "M_recoil is computed and stored in the ROOT ntuple for offline fitting; " \
  "a selection window is applied around the anti-Lambda_c- mass.")

alg.note(:deltaM_variable,
  "Delta_M = M(Lambda_c+ pi+ pi-) - M(Lambda_c+), where M(Lambda_c+) is " \
  "the invariant mass of the tagged Lambda_c+ candidate from the tag side. " \
  "Delta_M is a key discriminant: it peaks at m_Lc(2625)+ - m_Lc+ " \
  "(~0.341 GeV/c^2) for signal events. The best pi+ pi- combination is " \
  "selected by minimum |Delta_M - (m_Lc(2625)+ - m_Lc+)| when multiple " \
  "combinations are possible. This selection is applied in the ROOT " \
  "analysis stage, not within the BOSS DSL framework.")

alg.note(:vertex_fit,
  "A Lambda_c+ pi+ pi- vertex fit is performed before the final selection. " \
  "S_bachelor tracks (the pi+ and pi- from Lambda_c(2625)+ decay) and " \
  "S_daughter tracks (the Lambda_c+ decay daughters from the tag side: " \
  "p, K, pi for pK-pi+ mode; p, pi+, pi- for pK_S0 mode; p, pi-, pi+ for " \
  "Lambda pi+ mode) are constrained to a common vertex. The vertex fit " \
  "chi^2 is used as a selection criterion. This vertex fit is implemented " \
  "in the BOSS analysis module and is not directly expressible in the " \
  "TagAnalysis DSL.")

alg.note(:m_lc_pipi_and_mrecoil_fit,
  "Signal yields are extracted via a simultaneous unbinned maximum " \
  "likelihood fit to the M(Lambda_c+ pi+ pi-) and M_recoil distributions " \
  "in ROOT. The M(Lambda_c+ pi+ pi-) fit model includes: " \
  "Lambda_c(2625)+ signal (Crystal Ball function), Lambda_c(2595)+ " \
  "background (Gaussian), and combinatorial background (2nd-order " \
  "Chebyshev polynomial). Fits are performed separately for each energy " \
  "point and each tag mode. The M_recoil distribution is used to constrain " \
  "the background shapes. Branching fraction is computed as: " \
  "BF = N_sig / (N_tag * epsilon_sig), where N_tag is the single-tag yield " \
  "determined from fits to the M_BC distribution of tagged Lambda_c+ " \
  "candidates.")

alg.note(:tag_yield_determination,
  "Tag-side Lambda_c+ candidates are selected by the M_BC (beam-constrained " \
  "mass) distribution. For the pK-pi+ mode: |M_BC - m_Lc+| < 5 MeV/c^2. " \
  "Analogous mass windows applied for pK_S0 and Lambda pi+ modes. " \
  "Single-tag yields N_tag are determined by fitting the M_BC distribution " \
  "with a double-Gaussian signal shape plus an ARGUS background function, " \
  "performed separately for each tag mode and each energy point.")

alg.note(:decay_card,
  "Decay card uses 'vpho' (virtual photon) as the e+e- initial state " \
  "for MC generation at sqrt(s) = 4.918 and 4.950 GeV, which is the " \
  "standard BESIII convention for continuum energies above open-charm " \
  "threshold. The Lambda_c(2625)+ -> Lambda_c+ pi+ pi- decay is modelled " \
  "with PHSP (phase space). Lambda_c+ -> p+ K- pi+ is used as a " \
  "representative tag-mode decay in the card. anti-Lambda_c- -> anti-p- K+ pi- " \
  "decay is arbitrary as it cancels in the efficiency determination. " \
  "The charge-conjugate decay chain (e+e- -> Lambda_c+ anti-Lambda_c(2625)-) " \
  "is generated via a separate decay card. Signal MC cross-section line " \
  "shapes from BESIII measurements of e+e- -> anti-Lambda_c- Lambda_c(2625)+ " \
  "are applied as event weights in the ROOT analysis stage.")

alg.note(:best_combination,
  "When multiple pi+ pi- combinations exist in an event, the best " \
  "combination is selected by minimizing |Delta_M - (m_Lc(2625)+ - m_Lc+)|. " \
  "This selection is performed in the ROOT analysis stage and is not " \
  "expressible in the TagAnalysis DSL signal_side block.")

### Apply and execute ###
alg.with_decay_card(decay_card).apply
alg.execute_on(all_data + all_incMC + exMC_signal + exMC_cc)