### Dataset description ###
data_3773 = DatasetManager.real_data.find("712_3773")   # psi(3770), 2.92 fb^-1
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay card: psi(3770) -> D+ D-, one side D -> K_L0 e nu, other side hadronic tag
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                          PHSP;
  Enddecay

  Decay D+
  1.0000 anti-K0 e+ nu_e                PHOTOS ISGW2;
  Enddecay

  Decay anti-K0
  1.0000 K_L0                           PHSP;
  Enddecay

  Decay D-
  0.2000 K+ pi- pi-                     PHSP;
  0.2000 K+ pi- pi- pi0                 PHSP;
  0.2000 K_S0 pi- pi0                   PHSP;
  0.2000 K_S0 pi- pi- pi+               PHSP;
  0.1000 K_S0 pi-                       PHSP;
  0.1000 K+ K- pi-                      PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                    PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_KL0_e_nu"
  config.related_dataset = data_3773
  config.events = 500000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Tag Analysis (D-tag single tag, signal side: e+ + missing K_L0) ###
alg = TagAnalysis.new("DpKL0enu")
alg.set_header(["DpKL0enuAlg/DpKL0enu.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# ST tag on D+/- (both charges scanned by omitting charm); 6 hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKsPi, :DptoKKPi
end

# Signal side: one positron (opposite sign of the tag) + one shower giving K_L0 direction,
# with K_L0 declared as missing (massive form, resolved via U_miss = 0)
alg.signal_side do |s|
  s.charged(ep: 1, at_least: false)
  s.photons 1
  s.min_photon_energy 0.1              # K_L0 shower energy > 0.1 GeV
  s.min_photon_angle 10.0
  s.missing :K_L0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:signal_electron_pid,
         "Electron PID uses combined EMC + dE/dx + TOF likelihoods with " \
         "L'(e)/[L'(e)+L'(pi)+L'(K)] > 0.8; bremsstrahlung recovery adds " \
         "unmatched showers within 5 deg of the electron; the electron must " \
         "have opposite charge to the ST D and no other unused charged track is allowed.")
   .note(:kL0_shower_direction,
         "The K_L0 direction is taken from an unused neutral shower; its " \
         "momentum magnitude is solved from four-momentum conservation and " \
         "U_miss = 0. If multiple K_L0 shower candidates exist, the most " \
         "energetic shower is chosen.")
   .note(:pi0_veto_on_KL0_shower,
         "Reject photons that form part of a gamma-gamma pair with " \
         "0.110 < M(gg) < 0.155 GeV to suppress pi0 photons feeding the K_L0 shower.")
   .note(:tag_deltaE_windows,
         "Mode-dependent Delta E windows for ST D candidates (about +/-3 sigma): " \
         "K pi pi (+/-30 MeV), K pi pi pi0 (-52,+39 MeV), Ks pi pi0 (-57,+40 MeV), " \
         "Ks pi pi pi (+/-34 MeV), Ks pi (+/-32 MeV), K K pi (+/-30 MeV).")
   .note(:tag_mBC_signal_region,
         "Tag mBC signal region: 1.86 < M_BC < 1.88 GeV; only one candidate " \
         "per mode kept (smallest |Delta E|).")
   .note(:pi0_selection,
         "pi0 candidates from photon pairs with 0.110 < M(gg) < 0.155 GeV; " \
         "1C mass-constrained kinematic fit to nominal pi0 mass, chi2 < 20.")
   .note(:kS0_selection,
         "K_S0 from oppositely charged track pairs with vertex fit; vertex within " \
         "20 cm along beam axis; |M(pi+pi-) - m(K_S0)| < 12 MeV; decay length > 2 sigma from IP.")
   .with_decay_card(decay_card_signal)
   .apply

alg.execute_on([data_3773, incMC_3773, exMC_signal])
