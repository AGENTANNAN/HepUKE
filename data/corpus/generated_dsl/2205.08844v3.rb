### Dataset description ###
# Real data — six energy points spanning 4.178-4.226 GeV (total ~6.32 fb^-1)
data_4178 = DatasetManager.real_data.find("703_4180")   # <Ecm> ~ 4178 MeV
data_4190 = DatasetManager.real_data.find("703_4190")   # <Ecm> ~ 4188.8 MeV
data_4200 = DatasetManager.real_data.find("703_4200")   # <Ecm> ~ 4198.9 MeV
data_4210 = DatasetManager.real_data.find("703_4210")   # <Ecm> ~ 4209.2 MeV
data_4220 = DatasetManager.real_data.find("703_4220")   # <Ecm> ~ 4218.7 MeV
data_4226 = DatasetManager.real_data.find("703_4230")   # <Ecm> ~ 4226.26 MeV

# Corresponding inclusive MC samples
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
# No exclusive-MC sample is generated for this DSL.

# Decay card: e+e- -> D_s*+ D_s- -> gamma D_s+ D_s-
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s+
    1.0000 K+ pi+ pi- PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Event selection (BOSS, tag-based) ###
alg = TagAnalysis.new("DsDsStarTag")
alg.set_header(["DsDsStarTagAlg/DsDsStarTag.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })   # nominal value; the fit uses the per-run measured beam energy
   .with_decay_card(decay_card)

# Tag 1 — Ds- (charm -1), ten hadronic tag modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKPi0, :DstoKKPi, :DstoKPi0Pi0, :DstoKKPiPi0, :DstoKKPiPiPi,
          :DstoKsK, :DstoPiPiPi, :DstoPiEta, :DstoPiEtaPrime, :DstoKsEta
  t.charm -1
end

# Tag 2 — Ds+ (charm +1) through the signal mode K+ pi+ pi-, best combination by smallest invariant mass
alg.tag_side(:Ds) do |t|
  t.modes :DstoKPiPi
  t.charm 1
  t.rank_by :inv
end

# Signal side — exactly one transition photon, isolated from charged tracks
alg.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
end

# 6C/7C kinematic fit: 4-momentum conservation + Ds mass constraints + D_s* mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag1, :gamma).constrain_to_nominal_mass_of(:"D_s*+")
  f.invariant_mass_of(:tag2, :gamma).constrain_to_nominal_mass_of(:"D_s*+")
  f.chi2_cut 200
end

# BOSS-side cut not expressible in the tag DSL: K_S0 veto
alg.note(:background_veto,
  "reject events with M(pi+ pi-) in [0.4676, 0.5276] GeV/c^2 to suppress K_S0 contamination")

alg.apply
alg.execute_on([data_4178, data_4190, data_4200, data_4210, data_4220, data_4226,
                incMC_4178, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4226])