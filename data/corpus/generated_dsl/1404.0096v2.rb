# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# J/ψ (3.097 GeV) real data, 2009 + 2012 rounds (BOSS 708); sample 708_3097
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### Decay cards (EvtGen format) ###

# Mode A — signal: J/ψ → γ η′, η′ → π+π−π+π−
decay_card_modeA_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Mode A — background: J/ψ → γ η′, η′ → π+π−η, η → γπ+π−
decay_card_modeA_bkg_eta = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.000 gamma pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Mode A — background: J/ψ → γ f2(1270), f2(1270) → π+π−π+π−
decay_card_modeA_bkg_f2 = <<~DECAYCARD
    Decay J/psi
    1.000 gamma f2(1270) PHSP;
    Enddecay
    Decay f2(1270)
    1.000 pi+ pi- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Mode B — signal: J/ψ → γ η′, η′ → π+π−π0π0, π0 → γγ
decay_card_modeB_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- pi0 pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode B — background: J/ψ → γ η′, η′ → π+π−η, η → π0π0π0
decay_card_modeB_bkg_eta = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    1.000 pi0 pi0 pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode B — background: J/ψ → γ f2(1270), f2(1270) → π+π−π+π−
decay_card_modeB_bkg_f2 = <<~DECAYCARD
    Decay J/psi
    1.000 gamma f2(1270) PHSP;
    Enddecay
    Decay f2(1270)
    1.000 pi+ pi- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive MC (500k events each) ###
exMC_modeA_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_etap_2pip2pim"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeA_signal
  c.cross_section   = :default
end

exMC_modeA_bkg_eta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_etap_pipim_eta_gam2pi"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeA_bkg_eta
  c.cross_section   = :default
end

exMC_modeA_bkg_f2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_f2_4pi_charged_channel"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeA_bkg_f2
  c.cross_section   = :default
end

exMC_modeB_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_etap_pipim2pi0"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeB_signal
  c.cross_section   = :default
end

exMC_modeB_bkg_eta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_etap_pipim_eta_3pi0"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeB_bkg_eta
  c.cross_section   = :default
end

exMC_modeB_bkg_f2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "exmc_jpsi_gamma_f2_4pi_neutral_channel"
  c.related_dataset = jpsi_data
  c.events          = 500_000
  c.decay_card      = decay_card_modeB_bkg_f2
  c.cross_section   = :default
end

### Event selection — Mode A: J/ψ → γ π+π−π+π− (via η′ → π+π−π+π−) ###
algA_name = "JpsiGammaEtapCharged"
alg_modeA = Algorithm.new(algA_name)
alg_modeA.set_header(["#{algA_name}Alg/#{algA_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_modeA = Selection.new
sel_modeA.select_track {
      cos_theta 0.93    # |cosθ| < 0.93
      Vz 10.0           # |Vz| < 10 cm
      Vr 1.0            # Vr < 1 cm
      nChrp "==2"       # exactly two positive charged tracks
      nChrn "==2"       # exactly two negative charged tracks
      nNet "==0"        # zero net charge
    }
    .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      angle_to_track 10.0       # ≥ 10° from the nearest charged track
      energyThreshold_b 0.025   # E > 25 MeV in the EMC barrel
      energyThreshold_e 0.050   # E > 50 MeV in the EMC endcap
      nGam ">=1"                # at least one photon (radiative γ)
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # π+ and π− (charge-conjugation shorthand)
      npip ">=1"      # at least one identified π+
      npim ">=1"      # at least one identified π−
    }
    .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks assumed pions
    .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {   # nominal 4C γπ+π−π+π−
      nominal
      vertex_fit([1, 2, 3, 4])   # vertex fit on the four charged tracks
      constrain_four_momentum    # 4C energy-momentum constraint
      chi2_cut 200               # loose cut; tight χ² < 35 applied in ROOT
    }
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) {  # competing 4C γγπ+π−π+π−
      constrain_four_momentum    # stores the competing χ² for a ROOT-level ordering veto
    }

alg_modeA.with_decay_card(decay_card_modeA_signal).apply(sel_modeA)

### Event selection — Mode B: J/ψ → γ π+π−π0π0 (via η′ → π+π−π0π0) ###
algB_name = "JpsiGammaEtapNeutral"
alg_modeB = Algorithm.new(algB_name)
alg_modeB.set_header(["#{algB_name}Alg/#{algB_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_modeB = Selection.new
sel_modeB.select_track {
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
      nChrp "==1"       # exactly one positive charged track
      nChrn "==1"       # exactly one negative charged track
      nNet "==0"        # zero net charge
    }
    .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      angle_to_track 10.0       # ≥ 10° from the nearest charged track
      energyThreshold_b 0.025   # E > 25 MeV in the EMC barrel
      energyThreshold_e 0.050   # E > 50 MeV in the EMC endcap
      nGam ">=5"                # at least five photons (2×π0 → 4γ plus the radiative γ)
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # π+ and π− (charge-conjugation shorthand)
      npip "==1"      # exactly one identified π+
      npim "==1"      # exactly one identified π−
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C mass-constrained π0 fits on photon pairs
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 50
      npi0 ">=2"      # at least two π0 candidates
    }
    .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) {   # nominal 6C γπ+π−π0π0
      nominal
      vertex_fit([1, 2])         # vertex fit on the two charged tracks
      constrain_four_momentum    # 4C energy-momentum constraint (6C together with the π0 masses)
      chi2_cut 200               # loose cut; tight χ² < 35 applied in ROOT
    }
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0, :pi0]) {  # competing 6C γγπ+π−π0π0
      constrain_four_momentum    # stores the competing χ² for a ROOT-level ordering veto
    }

alg_modeB.with_decay_card(decay_card_modeB_signal).apply(sel_modeB)

### Execute on datasets ###
root_files_A = alg_modeA.execute_on([jpsi_data, jpsi_incMC, exMC_modeA_signal, exMC_modeA_bkg_eta, exMC_modeA_bkg_f2])
root_files_B = alg_modeB.execute_on([jpsi_data, jpsi_incMC, exMC_modeB_signal, exMC_modeB_bkg_eta, exMC_modeB_bkg_f2])