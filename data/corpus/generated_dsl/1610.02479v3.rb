# ==========================================================================
# ψ(3686) → γ χ_c1 ,  χ_c1 → η π+π−   (three independent η decay modes)
#   Mode I   : η → γγ          (visible final state: γ γ γ π+π−)
#   Mode II  : η → π+π−π0      (visible final state: γ γ γ π+π−π+π−)
#   Mode III : η → π0π0π0      (visible final state: γ γ γ γ γ γ γ π+π−)
# Independent final states + different selection -> one Algorithm each (Rule T1).
# Post-fit (ROOT-level) tighter chi2 cuts, eta mass window, radiative-photon
# energy window and sideband vetoes are NOT expressed here (out of BOSS scope).
# ==========================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(2S) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # generic ψ(2S) inclusive MC

# ---------------------------------------------------------------- decay cards
# ---- Mode I : η → γγ ----
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Mode II : η → π+π−π0 , π0 → γγ ----
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Mode III : η → π0π0π0 , π0 → γγ ----
decay_card_modeIII = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ------------------------------------------------------ exclusive signal MC
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_eta_gg"
  config.related_dataset = psip_data
  config.events          = 800_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_eta_pipimpi0"
  config.related_dataset = psip_data
  config.events          = 460_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_eta_3pi0"
  config.related_dataset = psip_data
  config.events          = 660_000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

### ======================= Event selection (BOSS) ======================= ###

# ----------------------------- Mode I : η → γγ -----------------------------
alg_name_I = "Chic1EtaModeI"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_I = Selection.new
sel_I.select_track {                 # charged tracks: 1 π+ , 1 π−, all treated as pions
        cos_theta 0.93
        Vz        20.0
        Vr        2.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      }
     .select_photon {                # ≥3 photons (radiative + 2 from η)
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"
      }
     .assign({:chrgp => :pip, :chrgn => :pim})   # all charged tracks assumed pions
     .kalman_kinematic_fit([:gamma, :gamma]) {   # γγ mass constraint to the η
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 15
        neta ">=1"
      }
     .kinematic_fit([:gamma, :eta, :pip, :pim]) { # nominal 5C: 4C + m(γγ)=m(η)
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_I.with_decay_card(decay_card_modeI).apply(sel_I)
alg_I.execute_on([psip_data, psip_incMC, exMC_modeI])

# -------------------------- Mode II : η → π+π−π0 ---------------------------
alg_name_II = "Chic1EtaModeII"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_II = Selection.new
sel_II.select_track {                # charged tracks: 2 π+ , 2 π−, all treated as pions
        cos_theta 0.93
        Vz        20.0
        Vr        2.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      }
      .select_photon {               # ≥3 photons (radiative + 2 from π0)
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"
      }
      .assign({:chrgp => :pip, :chrgn => :pim})   # all charged tracks assumed pions
      .kalman_kinematic_fit([:gamma, :gamma]) {   # γγ mass constraint to the π0
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      .kinematic_fit([:gamma, :pi0, :pip, :pip, :pim, :pim]) { # nominal 5C: 4C + m(γγ)=m(π0)
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_II.with_decay_card(decay_card_modeII).apply(sel_II)
alg_II.execute_on([psip_data, psip_incMC, exMC_modeII])

# ------------------------- Mode III : η → π0π0π0 ---------------------------
alg_name_III = "Chic1EtaModeIII"
alg_III = Algorithm.new(alg_name_III)
alg_III.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_III = Selection.new
sel_III.select_track {               # charged tracks: 1 π+ , 1 π−, all treated as pions
        cos_theta 0.93
        Vz        20.0
        Vr        2.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      }
       .select_photon {              # ≥7 photons (radiative + 6 from three π0)
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=7"
      }
       .assign({:chrgp => :pip, :chrgn => :pim})  # all charged tracks assumed pions
       .kalman_kinematic_fit([:gamma, :gamma]) {  # γγ mass constraints building the π0's
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=3"
      }
       .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pip, :pim]) { # nominal 7C: 4C + 3×m(γγ)=m(π0)
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_III.with_decay_card(decay_card_modeIII).apply(sel_III)
alg_III.execute_on([psip_data, psip_incMC, exMC_modeIII])