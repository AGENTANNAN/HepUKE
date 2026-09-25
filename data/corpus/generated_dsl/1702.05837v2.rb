### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773") # ψ(3770) inclusive MC

# Decay card — signal: ψ(3770) → D+ D−, D− → K+ π− π−, D+ → γ e+ νe
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 gamma e+ nu_e PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — dominant background: ψ(3770) → D+ D−, D− → K+ π− π−, D+ → π0 e+ νe, π0 → γ γ
decay_card_bkg = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 pi0 e+ nu_e PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal D+ → γ e+ νe
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DpToGammaENu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# 100k-event exclusive MC for the dominant D+ → π0 e+ νe background
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DpToPi0ENu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) — double tag: hadronic tag D− + signal D+ → γ e+ νe ###
alg_name = "DpToGammaENu"
my_alg   = TagAnalysis.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card_signal)

# Tag side: hadronic D− reconstructed in the six single-tag modes (charge-conjugate modes implied)
my_alg.tag_side(:Dminus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                               # pin the tagged D to D−
  t.window :mBC, min: 1.8628, max: 1.8788  # single-tag M_BC window (explicitly requested)
end

# Signal side: D+ → γ e+ νe — one electron, one radiative photon, and a massless missing neutrino
my_alg.signal_side do |s|
  s.charged(ep: 1)          # exactly one leftover track, identified as the electron
  s.require_charge 1        # total signal-side charge +1 (electron opposite the D− tag)
  s.photons 1               # at least one remaining photon
  s.min_photon_angle 10.0   # photon more than 10 degrees from the track
  s.min_photon_energy 0.010 # photon energy > 10 MeV
  s.missing :nu_e           # massless missing neutrino
end

# 4C kinematic fit of the tag, γ, e+ and νe to the lab energy, χ² < 200
my_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that cannot be expressed with the tag DSL
my_alg
  .note(:tag_track_selection, "tag D- charged tracks required with |cos(theta)| < 0.93, Vr < 1 cm, |Vz| < 10 cm")
  .note(:tag_pid, "pi/K separation in the tag by PID likelihoods: L(pi) > L(K) for pions, L(K) > L(pi) for kaons")
  .note(:ks0_selection, "K_S0 -> pi+pi- tag candidates: |cos(theta)| < 0.93, |Vz| < 20 cm, no Vr and no PID requirement, 0.487 < M(pi+pi-) < 0.511 GeV/c^2, decay length > 2 sigma")
  .note(:photon_selection, "tag photons: E > 25 MeV in the barrel (|cos(theta)| < 0.80) or E > 50 MeV in the endcap (0.84 < |cos(theta)| < 0.92), EMC shower time within 700 ns")
  .note(:pi0_selection, "pi0 -> gamma gamma tag candidates: 0.115 < M(gamma gamma) < 0.150 GeV/c^2 and a 1C mass-constrained fit; rejected if both photons are in the endcaps")
  .note(:tag_ranking, "when several tag candidates occur in the same mode, keep the one with the smallest |deltaE|")
  .note(:tag_deltae_windows, "mode-dependent single-tag |deltaE| windows in MeV (applied at the analysis level): K pi pi [-27,25]; K pi pi pi0 [-62,34]; Ks pi [-25,25]; Ks pi pi0 [-73,41]; Ks pi pi pi [-33,30]; K K pi [-23,20]")
  .note(:signal_electron_pid, "signal-side electron (opposite the tag) identified with dE/dx, TOF and EMC: L(e) > 0 and L(e)/(L(e)+L(pi)+L(K)) > 0.8")
  .note(:radiative_photon_selection, "the highest-energy remaining photon is taken as the radiative candidate; FSR photons with E > 50 MeV within a 5 degree cone around the electron are added; the radiative photon lateral moment is required to lie in (0.0, 0.3)")
  .note(:background_veto, "any photon pair forming a pi0 with a 1C-fit chi^2 < 20 is vetoed")
  # The fitted missing mass, U_miss and q^2 are stored automatically by the fit diagnostics.

my_alg.apply
root_files = my_alg.execute_on([data_3773, incMC_3773, exMC_signal, exMC_bkg])