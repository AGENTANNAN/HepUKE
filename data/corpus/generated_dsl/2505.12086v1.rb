### Dataset description ###
# ψ(3770), √s = 3.773 GeV
data_3773  = DatasetManager.real_data.find("712_3773")       # Real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # Corresponding inclusive MC sample

# Decay card for the signal process ψ(3770) → D+D−, D+ → π+ηη, D− → K+π−π−, η → γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ eta eta PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal channel
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DptoPiEtaEta"
  config.related_dataset = data_3773        # associate with the 3.773 GeV real dataset
  config.events          = 100_000          # 100k events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — tag-based analysis ###
alg_name = "DpToPiEtaEtaTag"
d_tag_analysis = TagAnalysis.new(alg_name)
d_tag_analysis.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({ "ECMS" => [:double, 3.773] })   # √s = 3.773 GeV
              .with_decay_card(decay_card_signal)
              .note(:background_veto, "residual combinatorial background suppressed
                by a BDTG classifier; training/applying the BDTG (input variables from
                the tag and signal side) is performed outside the BOSS selection")

# Tag side: the D− (charm −1), single tag, six hadronic modes.
# Store-not-cut: tag M_BC / ΔE are stored unconditionally; only the explicitly
# requested ST windows below are emitted as guards.
d_tag_analysis.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # D− → K+π−π−
          :DptoKPiPiPi0,     # D− → K+π−π−π0
          :DptoKsPi,         # D− → K_S π−
          :DptoKsPiPi0,      # D− → K_S π−π0
          :DptoKsPiPiPi,     # D− → K_S π−π−π+
          :DptoKKPi          # D− → K+K−π−
  t.charm -1                 # pin the tagged side to charm −1 (the D−)
  t.window :mBC,   min: 1.860, max: 1.880   # ST window M_BC ∈ [1.860, 1.880] GeV
  t.window :deltaE, abs: 0.040              # ST window |ΔE| < 0.040 GeV
end

# Signal side: the D+ → π+ηη (η → γγ) recoil against the tag.
# Exactly one π+ and four photons; minimum photon angle to charged tracks 10°.
d_tag_analysis.signal_side do |s|
  s.charged(pip: 1)         # exactly one π+
  s.photons 4               # four photons (two η → γγ)
  s.min_photon_angle 10.0   # photon–track opening angle > 10°
end

# 4C kinematic fit: energy-momentum conservation with both η masses constrained;
# χ² < 200. (The M(π+ηη) → nominal D+ mass constraint is applied at the ROOT stage,
# not here.)
d_tag_analysis.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # first η
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # second η
  f.chi2_cut 200
end

# Render the tag analysis (no Selection argument for TagAnalysis)
d_tag_analysis.apply

# Run the same selection over data, inclusive MC, and the signal exclusive MC
root_files = d_tag_analysis.execute_on([data_3773, incMC_3773, exMC_signal])