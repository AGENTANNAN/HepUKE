# frozen_string_literal: true
# ===========================================================================
# Semileptonic Ds+ branching fractions via e+e- -> Ds*+ Ds*-
# (hadronic tag Ds*-  +  signal-side Ds+ -> (eta/eta'/phi/f0/K0/K*0) e+ nu_e)
# BOSS part only: dataset preparation + tag-based (double-tag) selection.
# ===========================================================================

### Dataset preparation ###
# The 4.237-4.699 GeV scan; sample name convention "<BOSS version>_<CMS energy in MeV>"
scan_names = %w[
  703_4237 703_4246 703_4260 703_4270 703_4280 703_4310 703_4360 703_4390
  703_4420 703_4470 703_4530 703_4575 703_4600
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
]
scan_data = scan_names.map { |n| DatasetManager.real_data.find(n) }

# Matching inclusive MC samples (available for the subset of points below)
scan_incMC = %w[
  703_4237 703_4246 703_4260 703_4270 703_4280 703_4360 703_4420 703_4600
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
].map { |n| DatasetManager.inclusive_mc.find(n) }

# Signal decay card of the production process
#   psi(4260) -> Ds*+ Ds*-,  Ds*+- -> gamma Ds+-,  Ds+- -> K+ K- pi+-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive signal MC for every scan point (same card, one MC per point)
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_dsst_dsst"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ===========================================================================
# Mode 1 : Ds+ -> eta e+ nu_e,  eta -> gamma gamma
# ===========================================================================
decay_card_eta = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 eta e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

alg_eta = TagAnalysis.new("DsTagSLeta")
alg_eta.set_header(["DsTagSLetaAlg/DsTagSLeta.h"])
       .set_constant({ "ECMS" => [:double, 4.260] })

# Tag side: single tag Ds*- (Ds*- -> gamma(pi0) Ds-), charm = -1
alg_eta.tag_side(:Ds) do |t|
  t.mode_group :hadronic       # 14 hadronic Ds- tag modes
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123   # explicit M_BC signal region (energy dependent)
end

# Signal side: tracks/showers not used by the tag; one electron + missing neutrino
alg_eta.signal_side do |s|
  s.charged(ep: 1)             # exactly one signal-side electron, no extra charged tracks
  s.photons 2                  # eta -> gamma gamma
  s.missing :nu_e              # missing (massless) neutrino
  s.min_photon_angle 10.0      # bremsstrahlung recovery: EMC showers within 10 deg
end

alg_eta.fit do |f|
  f.constrain_four_momentum                       # 4C constraint against the measured CMS
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma, :ep).between(0.0, 1.75)   # Cabibbo-suppressed hadron-lepton mass
  f.chi2_cut 200
end

alg_eta.note(:pid_correction_method,
             "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8; " \
             "this CL-based PID is not DSL-tunable (tag-DSL lepton keys use fixed v1 thresholds)")
     .note(:background_veto,
           "double-tag side vetoed for extra energy: maximum unused EMC shower energy < 0.3 GeV; " \
           "no extra charged tracks beyond the declared signal/tag content")
     .note(:efficiency_curve,
           "M_BC signal region is energy dependent (2.104-2.123 GeV/c^2 at the Ds mass point); " \
           "ST yields are extracted from fits to the stored M_BC distribution and " \
           "M^2_miss fits use the Ds*+ momentum constrained to the tag Ds*- direction")

alg_eta.with_decay_card(decay_card_eta).apply
alg_eta.execute_on(scan_data + scan_incMC + exMC_signal)

# ===========================================================================
# Mode 2 : Ds+ -> eta' e+ nu_e,  eta' -> eta pi+ pi-,  eta -> gamma gamma
# ===========================================================================
decay_card_etap = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 eta' e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

alg_etap = TagAnalysis.new("DsTagSLetap")
alg_etap.set_header(["DsTagSLetapAlg/DsTagSLetap.h"])
        .set_constant({ "ECMS" => [:double, 4.260] })

alg_etap.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123
end

alg_etap.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.photons 2
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_etap.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:etap)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma, :ep).between(0.0, 1.75)  # eta and eta' fitted simultaneously
  f.chi2_cut 200
end

alg_etap.note(:pid_correction_method,
              "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8 " \
              "(CL-based PID not expressible in the tag DSL)")
       .note(:background_veto,
             "maximum unused EMC shower energy < 0.3 GeV and no extra charged tracks on the double-tag side")

alg_etap.with_decay_card(decay_card_etap).apply
alg_etap.execute_on(scan_data + scan_incMC + exMC_signal)

# ===========================================================================
# Mode 3 : Ds+ -> phi e+ nu_e,  phi -> K+ K-
# ===========================================================================
decay_card_phi = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 phi e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

alg_phi = TagAnalysis.new("DsTagSLphi")
alg_phi.set_header(["DsTagSLphiAlg/DsTagSLphi.h"])
       .set_constant({ "ECMS" => [:double, 4.260] })

alg_phi.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123
end

alg_phi.signal_side do |s|
  s.charged(kp: 1, km: 1, ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_phi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :km).between(1.004, 1.034)          # phi mass window
  f.invariant_mass_of(:kp, :km, :ep).between(0.0, 1.90)       # Cabibbo-favored hadron-lepton mass
  f.chi2_cut 200
end

alg_phi.note(:pid_correction_method,
             "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8 " \
             "(CL-based PID not expressible in the tag DSL)")
      .note(:background_veto,
            "maximum unused EMC shower energy < 0.3 GeV and no extra charged tracks on the double-tag side")

alg_phi.with_decay_card(decay_card_phi).apply
alg_phi.execute_on(scan_data + scan_incMC + exMC_signal)

# ===========================================================================
# Mode 4 : Ds+ -> f0(980) e+ nu_e,  f0(980) -> pi+ pi-
# ===========================================================================
decay_card_f0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 f_0 e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay f_0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

alg_f0 = TagAnalysis.new("DsTagSLf0")
alg_f0.set_header(["DsTagSLf0Alg/DsTagSLf0.h"])
      .set_constant({ "ECMS" => [:double, 4.260] })

alg_f0.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123
end

alg_f0.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_f0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.880, 1.080)        # f0(980) mass window
  f.invariant_mass_of(:pip, :pim, :ep).between(0.0, 1.90)
  f.chi2_cut 200
end

alg_f0.note(:pid_correction_method,
            "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8 " \
            "(CL-based PID not expressible in the tag DSL)")
     .note(:background_veto,
           "maximum unused EMC shower energy < 0.3 GeV and no extra charged tracks on the double-tag side")

alg_f0.with_decay_card(decay_card_f0).apply
alg_f0.execute_on(scan_data + scan_incMC + exMC_signal)

# ===========================================================================
# Mode 5 : Ds+ -> K0 e+ nu_e,  K0 -> K_S0 -> pi+ pi-
# ===========================================================================
decay_card_k0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 K_S0 e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

alg_k0 = TagAnalysis.new("DsTagSLK0")
alg_k0.set_header(["DsTagSLK0Alg/DsTagSLK0.h"])
      .set_constant({ "ECMS" => [:double, 4.260] })

alg_k0.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123
end

alg_k0.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_k0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim, :ep).between(0.0, 1.90)
  f.chi2_cut 200
end

alg_k0.note(:pid_correction_method,
            "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8 " \
            "(CL-based PID not expressible in the tag DSL)")
     .note(:background_veto,
           "maximum unused EMC shower energy < 0.3 GeV and no extra charged tracks on the double-tag side")
     .note(:efficiency_curve,
           "K0 -> K_S0 -> pi+pi- requires a displaced secondary vertex; no secondary-vertex primitive " \
           "is available on the tag signal side, the K_S0 flight length / vertex quality is handled downstream")

alg_k0.with_decay_card(decay_card_k0).apply
alg_k0.execute_on(scan_data + scan_incMC + exMC_signal)

# ===========================================================================
# Mode 6 : Ds+ -> K*0 e+ nu_e,  K*0 -> K+ pi-
# ===========================================================================
decay_card_kstar = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 K*0 e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay K*0
  1.000 K+ pi- VSS;
  Enddecay

  End
DECAYCARD

alg_kstar = TagAnalysis.new("DsTagSLKstar0")
alg_kstar.set_header(["DsTagSLKstar0Alg/DsTagSLKstar0.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })

alg_kstar.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm -1
  t.window :mBC, min: 2.104, max: 2.123
end

alg_kstar.signal_side do |s|
  s.charged(kp: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_kstar.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :pim).between(0.882, 0.992)        # K*0 mass window
  f.invariant_mass_of(:kp, :pim, :ep).between(0.0, 1.90)
  f.chi2_cut 200
end

alg_kstar.note(:pid_correction_method,
               "signal-side electron identified by CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8 " \
               "(CL-based PID not expressible in the tag DSL)")
        .note(:background_veto,
              "maximum unused EMC shower energy < 0.3 GeV and no extra charged tracks on the double-tag side")
        .note(:efficiency_curve,
              "the best tag candidate is chosen per mode by minimum |DeltaE|; DeltaE and M_BC are stored " \
              "unconditionally and the energy-dependent M_BC signal region is applied downstream")

alg_kstar.with_decay_card(decay_card_kstar).apply
alg_kstar.execute_on(scan_data + scan_incMC + exMC_signal)