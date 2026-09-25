### Dataset preparation ###
# ψ(3770) real data (7.93 fb^-1) and the corresponding inclusive MC (BOSS 712, √s = 3.773 GeV)
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process: ψ(3770) → D0 anti-D0,
# D0 → a0(980)- e+ nu_e, a0(980)- → eta pi-, eta → gamma gamma,
# tag side anti-D0 → K+ pi- (phase space).
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 a_0(980)- e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay a_0(980)-
    1.000 eta pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200,000 exclusive-MC events for the signal decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toa0enu_a0toetapi_etatogg"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (tag-based, BOSS) ###
alg_name = "D0A0EtaNu"
tag_analysis = TagAnalysis.new(alg_name)
tag_analysis.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            # Signal-side electron PID uses the fixed v1 tag-framework thresholds (SimplePIDSvc);
            # the intended E/pc > 0.8 cut and the combined e/pi/K likelihood criteria are not
            # expressible inside the tag framework and are handled outside BOSS.
            .note(:pid_correction_method, "Signal-side electron identification relies on the fixed v1
              tag-framework thresholds; the intended E/pc > 0.8 requirement and the combined
              electron/pion/kaon likelihood criteria cannot be encoded in the tag framework.")
            .with_decay_card(decay_card_signal)

# Tag side: anti-D0 reconstructed from pre-stored DTag candidates (single tag).
# Modes K pi, K pi pi0, K pi pi pi; charm = -1 pins the tagged (anti-D0) side.
tag_analysis.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035   # explicit tag-side ΔE window
end

# Signal side: what the tag did not use — two photons, exactly one pi- and one e+,
# net charge zero, and the undetected (massless) neutrino.
tag_analysis.signal_side do |s|
  s.photons 2
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# Kinematic fit: 4-momentum conservation plus the eta mass constraint on the two
# signal photons, with chi2 < 200.
tag_analysis.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

# Validate and render the tag specification (no Selection argument).
tag_analysis.apply

# Execute on real data, inclusive MC, and the signal exclusive MC.
root_files = tag_analysis.execute_on([data_3773, incMC_3773, exMC_signal])