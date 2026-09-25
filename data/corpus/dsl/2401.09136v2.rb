# ============================================================
# BESIII Analysis: eta/eta' -> gamma e+ e- Dalitz decays
# from J/psi -> gamma eta/eta' with 10B J/psi events
# Paper: 2401.09136v2 — Improved measurements of eta/eta' Dalitz decays
# Form factor measurement (single-pole for eta, multipole for eta')
# Dark photon search A' -> e+ e-
# ============================================================

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for eta -> gamma e+ e- signal
decay_card_eta = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for eta' -> gamma e+ e- signal
decay_card_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for eta -> gamma e+ e- channel
exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "eta_gamma_ee"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_eta
  config.cross_section = :default
end

# Exclusive MC for eta' -> gamma e+ e- channel
exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "etap_gamma_ee"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_etap
  config.cross_section = :default
end

### Event selection (BOSS) — shared by both eta and eta' channels ###
# Both channels have identical final state (gamma, gamma, e+, e-) and selection criteria
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz 10.0
                  Vr 1.0
                  nChrp "==1"
                  nChrn "==1"
                  nNet "==0"
                }
               .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end 14
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  angle_to_track 10.0
                  nGam ">=2"
                }
               .pid(method: :probability) {
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"; nlm "==1"
                }
               .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 100
                }

# Notes for inexpressible BOSS-side procedures
boss_notes = lambda do |alg|
  alg.note(:pid_correction_method, "electron PID: L(e)/(L(e)+L(pi)) > 0.5 for eta channel, > 0.95 for eta' channel. PID efficiency corrected with radiative Bhabha and J/psi->e+e- control samples; 2D correction as function of cos(theta) and transverse momentum")
     .note(:helix_correction, "track helix parameter correction applied to MC to match data; efficiency difference with/without correction taken as systematic uncertainty for 4C kinematic fit")
     .note(:background_veto, "post-4C-fit vetoes: photon conversion rejection (R_xy > 2 cm when cos_theta_eg > 0 and |Delta_xy| < 0.8 cm); M(gamma gamma) veto |M-0.547|>0.03 GeV and |M-0.958|>0.03 GeV to suppress J/psi->e+e- eta/eta' with eta/eta'->gamma gamma; low-E photon energy > 0.15 GeV; min angle between radiative photon and e+/e- > 20 deg; min angle between low-E photon and e+/e- > 10 deg")
     .note(:efficiency_curve, "detection efficiency as function of M(e+e-) used for TFF form factor measurement; signal MC generated with custom generator incorporating theoretical TFF amplitudes; form factor parameters determined iteratively")
end

### Algorithm for eta -> gamma e+ e- channel ###
alg_eta = Algorithm.new("EtaDalitz")
alg_eta.set_header(["EtaDalitzAlg/EtaDalitz.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
boss_notes.call(alg_eta)
alg_eta.with_decay_card(decay_card_eta).apply(event_selection)
alg_eta.execute_on([jpsi_data, jpsi_incMC, exMC_eta])

### Algorithm for eta' -> gamma e+ e- channel ###
alg_etap = Algorithm.new("EtapDalitz")
alg_etap.set_header(["EtapDalitzAlg/EtapDalitz.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
boss_notes.call(alg_etap)
alg_etap.with_decay_card(decay_card_etap).apply(event_selection)
alg_etap.execute_on([jpsi_data, jpsi_incMC, exMC_etap])