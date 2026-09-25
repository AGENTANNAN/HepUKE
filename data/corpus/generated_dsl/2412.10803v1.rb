# ============================================================================
#  psi(3770) -> D0 D0bar
#  single tag : anti-D0 in three hadronic modes
#  signal     : D0 -> Kbar0 pi- e+ nu_e,  Kbar0 -> K_S0 -> pi+ pi-
# ============================================================================

### ---------------------------- Datasets --------------------------------- ###
data_3773  = DatasetManager.real_data.find("712_3773")     # 7.9 fb-1 psi(3770) data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # inclusive MC at 3.773 GeV

### ------------------- Decay card for the signal process ----------------- ###
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 anti-K0 pi- e+ nu_e PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### -------------------- Exclusive signal MC : 500k events ---------------- ###
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_d0_kspimenu"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### ------------------------- Tag analysis (BOSS) ------------------------- ###
alg = TagAnalysis.new("D0TagKsPiEnu")
alg.set_header(["D0TagKsPiEnuAlg/D0TagKsPiEnu.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # c.m. energy
   .with_decay_card(decay_card_signal)

# ---- Tag side: single tag of the anti-D0 (three hadronic modes) -----------
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0   # anti-D0 -> K+pi-, K+pi-pi-pi+, K+pi-pi0
  t.charm -1                                    # pin the tagged side to the anti-D0
end

# ---- Signal side: D0 -> Kbar0(->K_S0->pi+pi-) pi- e+ nu_e ---------------
alg.signal_side do |s|
  s.photons 0                      # no photon participates in the fit
  s.min_photon_energy 0.25         # shower energy floor 0.25 GeV (pi0 suppression threshold)
  s.charged(pip: 1, pim: 2, ep: 1) # K_S0 -> pi+pi- (1pi+ ,1pi-) , bachelor pi- , electron
  s.require_charge 0               # net charge zero
  s.missing :nu_e                  # one missing (massless) neutrino
end

# ---- Kinematic fit: 4C constraint over tag + signal + missing neutrino ----
alg.fit do |f|
  f.constrain_four_momentum                             # tag + signal + nu_e = c.m. four-momentum
  f.invariant_mass_of(:pip, :pim).between(0.485, 0.510) # K_S0 mass window (GeV/c^2)
  f.chi2_cut 200                                        # loose cut; tightened at ROOT level
end

### ---------- BOSS-side procedures not expressible in the DSL -------------- ###
alg.note(:ks0_reconstruction,
         "K_S0 reconstructed from the signal-side pi+pi- pair with decay-length "
       + "significance > 2 sigma; the secondary-vertex / decay-length criterion is "
       + "not expressible in the tag DSL (only the 0.485-0.510 GeV/c^2 mass window "
       + "is imposed in the fit).")
   .note(:electron_pid_selection,
         "Electron identified by L'_e > 0.001 and L'_e/(L'_e+L'_pi+L'_K) > 0.8, "
       + "together with E/p > 0.7 and bremsstrahlung recovery within 5 deg; the "
       + "signal-side ep key uses SimplePIDSvc with fixed thresholds and does not "
       + "encode these criteria.")
   .note(:background_veto,
         "Events with a photon of E_gamma > 0.25 GeV are vetoed to suppress the "
       + "pi0 background; the wrong-sign decay D0 -> Kbar0 pi+pi- is vetoed by "
       + "requiring M(Kbar0 pi- e+) < 1.80 GeV/c^2.")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])