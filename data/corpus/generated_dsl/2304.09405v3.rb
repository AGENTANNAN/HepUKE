# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Real data and inclusive MC at 4.600 GeV (reference point of the 7-point scan
# 4.600 / 4.612 / 4.628 / 4.641 / 4.661 / 4.682 / 4.699 GeV).
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card — reference mode: Lambda_c+ -> Sigma+ pi+ pi-, Sigma+ -> p pi0, pi0 -> gamma gamma
decay_card_ref = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ pi+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card — signal mode: Lambda_c+ -> Sigma+ K+ K-
decay_card_kk = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ K+ K- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive MC for the reference mode Sigma+ pi+ pi-
exMC_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_lambdac_ref_sigmapippi"
  config.related_dataset = data_4600
  config.events          = 500_000
  config.decay_card      = decay_card_ref
  config.cross_section   = :default
end

# 500k-event exclusive MC for the signal mode Sigma+ K+ K-
exMC_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_lambdac_sigmakk"
  config.related_dataset = data_4600
  config.events          = 500_000
  config.decay_card      = decay_card_kk
  config.cross_section   = :default
end

### Event selection (tag-based: single-tag Lambda_c+ from psi(4260) -> Lambda_c+ Lambda_c-) ###
alg_name = "LambdacST"
lam_tag = TagAnalysis.new(alg_name)                 # tag-based analysis layer (DTagTool driven)
lam_tag.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})

# Tag side: single tag on the Lambda_c+, reconstructed in the reference and signal modes
lam_tag.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoSigmaPiPi,     # Sigma+ pi+ pi-   (reference)
          :LambdacPtoSigmaKK,       # Sigma+ K+ K-
          :LambdacPtoSigmaKPi,      # Sigma+ K+ pi-
          :LambdacPtoSigmaKPiPi0,   # Sigma+ K+ pi- pi0
          :LambdacPtoSigmaPhi       # Sigma+ phi (phi -> K+ K-)
  t.charm 1                         # pin the Lambda_c+ (tag) charge
end

# Signal side: exactly two photons (the pi0 daughters), one K+, one K-, one proton, net charge +1
lam_tag.signal_side do |s|
  s.photons 2
  s.charged(kp: 1, km: 1, prp: 1)
  s.require_charge(1)
  s.min_photon_angle 10.0
end

# Kinematic fit: 4C energy-momentum conservation + gamma-gamma constrained to the nominal pi0 mass
lam_tag.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :gamma, :gamma).between(1.174, 1.200)   # Sigma+ : M(p pi0) window
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)         # pi0    : M(gamma gamma) window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Lambda / K_S0 vetoes have no DSL form on a tag fit — captured for the uncertainty layer
lam_tag.note(:background_veto,
  "veto Lambda: M(p pi-) outside [1.11, 1.12] GeV/c^2; veto K_S0: M(pi+ pi-) outside [0.48, 0.52] GeV/c^2")

# Mode-dependent DeltaE windows and the unbinned ML fits to M_BC are applied in the ROOT
# stage on the stored DeltaE / M_BC distributions (store-not-cut).

lam_tag.with_decay_card(decay_card_ref).apply   # tag apply takes NO Selection argument
root_files = lam_tag.execute_on([data_4600, incMC_4600, exMC_ref, exMC_kk])