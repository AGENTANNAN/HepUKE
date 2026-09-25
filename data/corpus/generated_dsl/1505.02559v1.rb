# =============================================================================
# Dataset preparation — untagged ISR analysis
#   e+e- -> gamma_ISR X(3872) -> gamma_ISR pi+pi- J/psi  (J/psi -> e+e- or mu+mu-)
#   the same final state also serves as the psi(3686) reference for Gamma_ee.
#   Four energy points: 4.009 / 4.230 / 4.260 / 4.360 GeV.
# =============================================================================
data_4009 = DatasetManager.real_data.find("703_4009")   # 4.009 GeV real data
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.230 GeV real data
data_4260 = DatasetManager.real_data.find("703_4260")   # 4.260 GeV real data
data_4360 = DatasetManager.real_data.find("703_4360")   # 4.360 GeV real data

incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")   # inclusive MC at 4.009 GeV
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")   # inclusive MC at 4.230 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")   # inclusive MC at 4.260 GeV
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")   # inclusive MC at 4.360 GeV

data_points   = [data_4009, data_4230, data_4260, data_4360]
inc_mc_points = [incMC_4009, incMC_4230, incMC_4260, incMC_4360]

# --- Decay cards (EvtGen syntax). psi(4260) is the KKMC top-mother convention ----
# Signal: ISR X(3872) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_x3872_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X(3872)  PHSP;
  Enddecay
  Decay X(3872)
  1.000  pi+ pi- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.000  e+ e-  PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Signal: ISR X(3872) -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_x3872_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma X(3872)  PHSP;
  Enddecay
  Decay X(3872)
  1.000  pi+ pi- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.000  mu+ mu-  PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Reference: ISR psi(2S) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_psi2S_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma psi(2S)  PHSP;
  Enddecay
  Decay psi(2S)
  1.000  pi+ pi- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.000  e+ e-  PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Reference: ISR psi(2S) -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_psi2S_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma psi(2S)  PHSP;
  Enddecay
  Decay psi(2S)
  1.000  pi+ pi- J/psi  PHSP;
  Enddecay
  Decay J/psi
  1.000  mu+ mu-  PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Background: e+e- -> eta J/psi, eta -> pi+ pi- pi0 (pi0 -> gamma gamma)
decay_card_bkg_etaJpsi = <<~DECAYCARD
  Decay psi(4260)
  1.000  eta J/psi  PHSP;
  Enddecay
  Decay eta
  1.000  pi+ pi- pi0  PHSP;
  Enddecay
  Decay pi0
  1.000  gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

# Background: e+e- -> gamma_ISR pi+ pi- pi+ pi-
decay_card_bkg_4pi = <<~DECAYCARD
  Decay psi(4260)
  1.000  gamma pi+ pi- pi+ pi-  PHSP;
  Enddecay
  End
DECAYCARD

# --- Exclusive MC (100k events each), sampled at every energy point ---------------
exMC_x3872_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_isr_x3872_ee"
  config.events        = 100000
  config.decay_card    = decay_card_x3872_ee
  config.cross_section = :default
end

exMC_x3872_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_isr_x3872_mumu"
  config.events        = 100000
  config.decay_card    = decay_card_x3872_mumu
  config.cross_section = :default
end

exMC_psi2S_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_isr_psi2S_ee"
  config.events        = 100000
  config.decay_card    = decay_card_psi2S_ee
  config.cross_section = :default
end

exMC_psi2S_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_isr_psi2S_mumu"
  config.events        = 100000
  config.decay_card    = decay_card_psi2S_mumu
  config.cross_section = :default
end

exMC_bkg_etaJpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_bkg_etaJpsi"
  config.events        = 100000
  config.decay_card    = decay_card_bkg_etaJpsi
  config.cross_section = :default
end

exMC_bkg_4pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_bkg_4pi_isr"
  config.events        = 100000
  config.decay_card    = decay_card_bkg_4pi
  config.cross_section = :default
end

# =============================================================================
# Event selection (BOSS)
# =============================================================================
# Signal algorithm (ISR X(3872)); the two lepton channels share one selection chain.
alg_x3872 = Algorithm.new("ISRX3872")
alg_x3872.set_header(["ISRX3872Alg/ISRX3872.h"])
         .set_constant({"ECMS" => [:double, 4.260]})         # nominal; beam energy is per-run
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_multienergy, "one algorithm is applied to four energy points " \
               "(4.009/4.230/4.260/4.360 GeV); the ECMS constant is only a nominal value and " \
               "the per-run beam energy must be taken from the framework")
         .note(:isr_photon_region, "the ISR photon is radiated nearly collinear with the beam " \
               "(|cos(theta_ISR)| > 0.95) and escapes undetected; its four-momentum is recovered " \
               "from the missing momentum in the 2C fit")

# Reference algorithm (ISR psi(3686)); same final state and selection as the signal.
alg_psi2S = Algorithm.new("ISRpsi2S")
alg_psi2S.set_header(["ISRpsi2SAlg/ISRpsi2S.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_multienergy, "one algorithm is applied to four energy points " \
               "(4.009/4.230/4.260/4.360 GeV); the ECMS constant is only a nominal value and " \
               "the per-run beam energy must be taken from the framework")
         .note(:isr_photon_region, "the ISR photon is radiated nearly collinear with the beam " \
               "(|cos(theta_ISR)| > 0.95) and escapes undetected; its four-momentum is recovered " \
               "from the missing momentum in the 2C fit")

# Common selection chain (shared by the e+e- and mu+mu- channels, and by the reference).
event_selection = Selection.new
    .select_track {                       # exactly four good charged tracks
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        1.0                     # Vr < 1 cm
        nChrp    "==2"                    # two positive tracks
        nChrn    "==2"                    # two negative tracks
        nNet     "==0"                    # net charge zero
    }
    .pid(method: :probability) {          # PID with the probability method
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.4  # p>1 GeV/c -> lepton; EMC E>0.4 GeV -> e, else mu
        identify :pion, against: [:kaon, :proton]   # pi+ and pi- (charge-conjugation shorthand)
        nlp  "==1"                        # one l+
        nlm  "==1"                        # one l-
        npip "==1"                        # one pi+
        npim "==1"                        # one pi-
    }
    .remove(:pip) { condition "three_momentum_of(:pip) > 0.6" }  # drop pion candidates above 600 MeV/c
    .remove(:pim) { condition "three_momentum_of(:pim) > 0.6" }
    .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {   # 2C fit on gamma_ISR pi+pi- l+l-
        nominal
        miss_track_of(:gamma)                          # missing ISR photon (massless) recovered from missing momentum
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # M(l+l-) = M(J/psi)
        chi2_cut 200                                   # loose BOSS-level cut (tight chi2 < 15 applied in ROOT)
    }

# Attach the decay cards and render the shared selection chain.
alg_x3872.with_decay_card(decay_card_x3872_ee).apply(event_selection)
alg_psi2S.with_decay_card(decay_card_psi2S_ee).apply(event_selection.dup)

# ==== Execute on real data, inclusive MC, signal/reference MC and background MC ====
all_datasets = data_points + inc_mc_points +
               exMC_x3872_ee + exMC_x3872_mumu +
               exMC_psi2S_ee + exMC_psi2S_mumu +
               exMC_bkg_etaJpsi + exMC_bkg_4pi

root_files_x3872 = alg_x3872.execute_on(all_datasets)
root_files_psi2S = alg_psi2S.execute_on(all_datasets)