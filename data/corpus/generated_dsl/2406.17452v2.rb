# frozen_string_literal: true
### Dataset preparation ###
# Seven real-data energy points, 4.128-4.226 GeV (705 @ 4.130/4.160; 703 @ 4.180-4.220)
data_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220")
]

# Matching inclusive MC samples (one per energy point)
incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220")
]

# Decay card for the signal process (EvtGen format, EvtGen particle names):
# e+e- -> Ds*+ Ds*-, Ds* -> gamma Ds (93.5%) or pi0 Ds (6.5%),
# signal Ds+ -> pi+ pi+ pi- pi0, pi0 -> gamma gamma.
# (The tag-side Ds- is left to decay generically through the default EvtGen table.)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Ds*+ Ds*- PHSP;
  Enddecay

  Decay Ds*+
  0.935 gamma Ds+ PHSP;
  0.065 pi0  Ds+ PHSP;
  Enddecay

  Decay Ds*-
  0.935 gamma Ds- PHSP;
  0.065 pi0  Ds- PHSP;
  Enddecay

  Decay Ds+
  1.000 pi+ pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC generated at each of the seven energy points
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsstar_dsstar"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — double-tag analysis ###
alg_name = "DsStarDsStarDTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })

# Inexpressible BOSS-side facts handed to the downstream skills
alg.note(:ecms_per_point, "the spec runs over seven energy points from 4.128 to 4.226 GeV; " \
                          "the ECMS constant is set to the dominant 4180 point (4.178 GeV) and " \
                          "the per-point beam energy is taken from each dataset at job time")
   .note(:alt_fit_8c, "an 8C variant of the kinematic fit that additionally constrains the Ds+ " \
                      "signal invariant mass is used for the amplitude analysis stage; the " \
                      "BOSS-level selection uses the 7C fit (4C + pi0 + tag Ds- + signal Ds+ mass) " \
                      "with chi2 < 200")

# Tag side: Ds- reconstructed in seven single-tag modes (charm -1 pins the Ds- member)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S K
          :DstoKKPi,       # K K pi
          :DstoKKPiPi0,    # K K pi pi0
          :DstoKsKPiPi,    # K_S K pi pi
          :DstoPiEta,      # pi eta (eta -> gamma gamma)
          :DstoPiEtaPrime, # pi eta' (eta' -> pi pi eta / gamma gamma)
          :DstoKPiPi       # K pi pi
  t.charm -1
end

# Signal side: Ds+ -> pi+ pi+ pi- pi0 (pi0 -> gamma gamma), net charge +1
alg.signal_side do |s|
  s.photons 2..2          # exactly two good showers from the pi0
  s.min_photon_angle 10.0 # photon opening angle with respect to charged tracks > 10 degrees
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
end

# 7C kinematic fit: 4-momentum conservation + pi0, tag-Ds- and signal-Ds+ mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:pip, :pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

# Tag-side M_BC > 1.85 GeV/c^2 and the D0 -> K- pi+, D+ -> K- pi+ pi+ vetoes are
# deliberately NOT applied here (store-not-cut): mBC / deltaE are stored and the
# windows are applied in the ROOT stage. No tag-side window is declared.
alg.apply

# Execute on all seven real-data points, their inclusive MC and the signal MC
root_files = alg.execute_on(data_points + incMC_points + exMCs)