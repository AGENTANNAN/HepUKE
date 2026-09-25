# Paper: 2406.17452v2
# Amplitude analysis of Ds+ → π+π+π-π0
# DT tag-based: tag Ds- from 7 modes → signal Ds+ → π+π+π-π0
# Energy range: 4.128-4.226 GeV (BOSS 705 + 703)
# 7.33 fb⁻¹ total integrated luminosity

# Datasets
ds_705_4130 = DatasetManager.real_data.find("705_4130")
ds_705_4160 = DatasetManager.real_data.find("705_4160")
ds_703_4180 = DatasetManager.real_data.find("703_4180")
ds_703_4190 = DatasetManager.real_data.find("703_4190")
ds_703_4200 = DatasetManager.real_data.find("703_4200")
ds_703_4210 = DatasetManager.real_data.find("703_4210")
ds_703_4220 = DatasetManager.real_data.find("703_4220")

incMC_705_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_705_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_703_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_703_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_703_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_703_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_703_4220 = DatasetManager.inclusive_mc.find("703_4220")

data_points = [ds_705_4130, ds_705_4160, ds_703_4180, ds_703_4190,
               ds_703_4200, ds_703_4210, ds_703_4220]
incMC_points = [incMC_705_4130, incMC_705_4160, incMC_703_4180, incMC_703_4190,
                incMC_703_4200, incMC_703_4210, incMC_703_4220]

# Decay card: e+e- → Ds*+ Ds*-, Ds*+ → γ/π0 Ds+, Ds+ → π+π+π-π0
# Tag side: Ds- reconstructed via DTagAlg; signal side: Ds+ → 2π+π-π0
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    0.935   gamma D_s+ PHSP;
    0.065   pi0   D_s+ PHSP;
    Enddecay

    Decay D_s*-
    0.935   gamma D_s- PHSP;
    0.065   pi0   D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0000  pi+ pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "Ds_to_3pip_pi0"
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# TagAnalysis for Ds+ → π+π+π-π0 via DT
# One tag_side for Ds- (7 ST modes), signal_side for Ds+ → 2π+π-π0
alg = TagAnalysis.new("DsTo3PiPi0_DT")
alg.set_header(["DsTo3PiPi0Alg/DsTo3PiPi0.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })
   .note(:ecms_per_dataset,
     "ECMS varies per dataset (4.128-4.226 GeV); per-run MeasuredEcmsSvc handles the actual value")

# Tag side: Ds- reconstruction via 7 ST tag modes
# Charm -1 pins the anti-Ds (Ds-) side
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoPiEtaGamGam, :DstoPiEtapPiPiEtaGamGam, :DstoKPiPi
  t.charm -1
end

# Signal side: Ds+ → π+π+π-π0 (2π+, 1π-, π0→γγ)
alg.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.charged(pip: 2, pim: 1)
  s.require_charge 1   # 2×(+) + 1×(-) = +1
end

# 7C kinematic fit: 4-momentum + π0 mass + Ds- tag mass + Ds+ signal mass
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:pip, :pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

alg
  .note(:amplitude_analysis,
    "Amplitude analysis with 11 intermediate components performed in ROOT;
     8C fit (adds Ds+ signal mass constraint) for amplitude analysis phase")
  .note(:ds_star_transition,
    "Ds*+ transition photon/π0 handled in ROOT analysis stage;
     Ds*+ mass constraint applied via Ds*+ candidate reconstruction")
  .note(:d0_dplus_veto,
    "D0→K-π+ and D+→K-π+π+ vetoes applied in ROOT to suppress backgrounds
     from incorrectly tagged D mesons")
  .note(:mBC_window,
    "Tag-side M_BC > 1.85 GeV/c² applied in ROOT stage")
  .with_decay_card(decay_card)
  .apply

alg.execute_on(data_points + incMC_points + exMCs)