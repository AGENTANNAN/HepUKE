### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # Corresponding inclusive MC sample

# Decay card for signal chain 1: D0 → ω γ′ (ω → π+π−π0, π0 → γγ); γ′ invisible, unconstrained mass
# The tag D̄0 decays via the three tagged modes K+π−, K+π−π0, K+π−π+π−
decay_card_omega_gammap = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 omega gamma_prime PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi- VSS;
    0.3333 K+ pi- pi0 PHSP;
    0.3333 K+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for signal chain 2: D0 → γ γ′ ; γ′ invisible, unconstrained mass
decay_card_gamma_gammap = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 gamma gamma_prime PHSP;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi- VSS;
    0.3333 K+ pi- pi0 PHSP;
    0.3333 K+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each of the two signal decay chains
exMC_omega_gammap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0_omega_gammap"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_omega_gammap
  config.cross_section   = :default
end

exMC_gamma_gammap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0_gamma_gammap"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_gamma_gammap
  config.cross_section   = :default
end

### Tag analysis — Mode A: D0 → ω γ′ ###
alg_name_omega = "D0ToOmegaGammap"
alg_omega = TagAnalysis.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})   # CMS energy = 3.773 GeV
         .with_decay_card(decay_card_omega_gammap)

# Tag side: single tag with D̄0 → K+π−, K+π−π0, K+π−π+π−
alg_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: ω γ′ (ω → π+π−π0 → π+π−γγ, γ′ invisible with unconstrained mass)
alg_omega.signal_side do |s|
  s.photons 2                                             # exactly two photons from the π0
  s.charged(pip: 1, pim: 1)                               # one π+ and one π−
  s.min_photon_angle 10.0                                 # photon angle > 10°
  s.min_photon_energy 0.025                               # photon energy > 25 MeV
  s.missing :gamma_prime, mass: nil                       # invisible γ′, unconstrained (massless form)
end

# Kinematic fit: four-momentum conservation + π0 and ω nominal-mass constraints, χ² < 200
alg_omega.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)              # π0 from γγ
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).between(0.700, 0.850)  # ω from π+π−γγ
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:omega)
  f.chi2_cut 200
end

alg_omega.note(:background_veto, "extra photons vetoed: total energy of all photons not used by the signal required below 0.1 GeV")
alg_omega.note(:endcap_veto,     "endcap background suppressed: |cos(theta_recoil)| < 0.7")

alg_omega.apply
alg_omega.execute_on([psi3770_data, psi3770_incMC, exMC_omega_gammap])

### Tag analysis — Mode B: D0 → γ γ′ ###
alg_name_gamma = "D0ToGammaGammap"
alg_gamma = TagAnalysis.new(alg_name_gamma)
alg_gamma.set_header(["#{alg_name_gamma}Alg/#{alg_name_gamma}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})   # CMS energy = 3.773 GeV
         .with_decay_card(decay_card_gamma_gammap)

# Tag side: single tag with D̄0 → K+π−, K+π−π0, K+π−π+π−
alg_gamma.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: exactly one photon plus the missing γ′ (invisible, unconstrained mass)
alg_gamma.signal_side do |s|
  s.photons 1                                             # exactly one photon
  s.min_photon_angle 10.0                                 # photon angle > 10°
  s.min_photon_energy 0.5                                 # photon energy > 0.5 GeV
  s.missing :gamma_prime, mass: nil                       # invisible γ′, unconstrained (massless form)
end

# Kinematic fit: four-momentum conservation, χ² < 200
alg_gamma.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_gamma.note(:background_veto, "extra photons vetoed: total energy of all photons not used by the signal required below 0.1 GeV")
alg_gamma.note(:endcap_veto,     "endcap background suppressed: |cos(theta_recoil)| < 0.7")

alg_gamma.apply
alg_gamma.execute_on([psi3770_data, psi3770_incMC, exMC_gamma_gammap])