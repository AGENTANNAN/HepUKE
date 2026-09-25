# BESIII measurement of D+ → μ+ ν_μ using DT method
# ArXiv: 2410.07626v2
# Dataset: 20.3 fb^-1 at √s = 3.773 GeV (ψ(3770))

### Dataset preparation ###
psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-  PHSP;
    Enddecay

    Decay D+
    1.0000  mu+  nu_mu  SLN;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_mu_nu"
  config.related_dataset = psip3770_data
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Algorithm: semileptonic ST tag — D+ → μ+ ν_μ ###
alg = TagAnalysis.new("DpMuNu")
alg.set_header(["DpMuNuAlg/DpMuNu.h"])
  .set_constant(ECMS: 3.773)
  .with_decay_card(decay_card)

# Tag side: ST D- reconstructed in 8 hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0,
          :DptoKsPiPiPi, :DptoKKPi, :DptoPiPiPi, :DptoKPiPiPiPi
  t.charm -1
end

# Signal side: D+ → μ+ ν_μ — one μ+, no photons, missing neutrino
alg.signal_side do |s|
  s.photons 0
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu
end

# Fit: tag + μ+ + ν = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:muon_pid, "muon PID uses hit depth d_μ+ vs p_μ+ and cosθ; EMC energy of muon candidate within (0.00, 0.35) GeV; no extra charged track allowed; maximum energy of unused shower E_max^{extra γ} < 0.3 GeV")
  .apply

alg.execute_on([psip3770_data, psip3770_incMC, exMC_signal])