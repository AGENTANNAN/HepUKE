# =====================================================================
# ψ(3770) single-tag measurement of absolute D-meson branching fractions
# BOSS part: dataset preparation + event selection up to the final 4C fit
# =====================================================================

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) data, √s = 3.773 GeV (2.93 fb^-1)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

### Decay cards (EvtGen) ###
# ψ(3770) → D Dbar; the tagged D decays to the signal mode, the opposite D to the
# charge-conjugate mode (charge-conjugate tags are naturally included).

# D0 → K- π+
decay_card_D0_KPi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# D0 → π+ π-
decay_card_D0_PiPi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# D+ → K+ π0, π0 → γγ
decay_card_Dp_KPi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K+ pi0 PHSP;
    Enddecay

    Decay D-
    1.0000 K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# D+ → K_S0 π+, K_S0 → π+ π-
decay_card_Dp_KsPi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K_S0 pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K_S0 pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# D0 → K_S0 π0, K_S0 → π+ π-, π0 → γγ
decay_card_D0_KsPi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K_S0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# D+ → π+ η, η → γγ
decay_card_Dp_PiEta = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ eta PHSP;
    Enddecay

    Decay D-
    1.0000 pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (100k events each) ###
exMC_D0_KPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_D0ToKPi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_D0_KPi
  config.cross_section   = :default
end

exMC_D0_PiPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_D0ToPiPi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_D0_PiPi
  config.cross_section   = :default
end

exMC_Dp_KPi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DpToKPi0"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_Dp_KPi0
  config.cross_section   = :default
end

exMC_Dp_KsPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DpToKsPi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_Dp_KsPi
  config.cross_section   = :default
end

exMC_D0_KsPi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_D0ToKsPi0"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_D0_KsPi0
  config.cross_section   = :default
end

exMC_Dp_PiEta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DpToPiEta"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_Dp_PiEta
  config.cross_section   = :default
end

# =====================================================================
### Event selection (BOSS) ###
# =====================================================================

# ---------------------------------------------------------------------
# Mode 1: D0 → K- π+   (two opposite-charge tracks, net charge zero)
# ---------------------------------------------------------------------
alg_name_D0_KPi = "D0ToKPiTag"
alg_D0_KPi = Algorithm.new(alg_name_D0_KPi)
alg_D0_KPi.set_header(["#{alg_name_D0_KPi}Alg/#{alg_name_D0_KPi}.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
          .note(:charge_conjugation, "D0 and anti-D0 tags are both accepted; the charge-conjugate " \
                "mode D0bar -> K+ pi- is covered by the same selection")

sel_D0_KPi = Selection.new
  .select_track {                    # charged track selection
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"                  # two opposite-charge tracks
    nChrn     "==1"
    nNet      "==0"                  # net charge zero
  }
  .pid(method: :probability) {       # PID: probability method, 0.001 cut
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K-
    identify :pion, against: [:kaon, :proton]   # pi+ and pi-
    nkm   ">=1"
    npip  ">=1"
  }
  .kinematic_fit([:km, :pip]) {      # final fit of the tagged D
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_D0_KPi.with_decay_card(decay_card_D0_KPi).apply(sel_D0_KPi)

# ---------------------------------------------------------------------
# Mode 2: D0 → π+ π-   (two opposite-charge tracks, net charge zero)
# ---------------------------------------------------------------------
alg_name_D0_PiPi = "D0ToPiPiTag"
alg_D0_PiPi = Algorithm.new(alg_name_D0_PiPi)
alg_D0_PiPi.set_header(["#{alg_name_D0_PiPi}Alg/#{alg_name_D0_PiPi}.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                 "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
           .note(:charge_conjugation, "D0 and anti-D0 tags both accepted")

sel_D0_PiPi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip  ">=1"
    npim  ">=1"
  }
  .kinematic_fit([:pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_D0_PiPi.with_decay_card(decay_card_D0_PiPi).apply(sel_D0_PiPi)

# ---------------------------------------------------------------------
# Mode 3: D+ → K+ π0,  π0 → γγ   (≥1 track, ≥2 photons)
# ---------------------------------------------------------------------
alg_name_Dp_KPi0 = "DpToKPi0Tag"
alg_Dp_KPi0 = Algorithm.new(alg_name_Dp_KPi0)
alg_Dp_KPi0.set_header(["#{alg_name_Dp_KPi0}Alg/#{alg_name_Dp_KPi0}.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                 "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
           .note(:pi0_mass_window, "reconstructed π0 required in the γγ mass window 0.115-0.150 GeV/c^2")

sel_Dp_KPi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"                  # at least one track
  }
  .select_photon {                   # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"          # at least two photons for π0 → γγ
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp  ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # π0 from γγ, mass constrained
    invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :pi0]) {      # final fit of the tagged D
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_Dp_KPi0.with_decay_card(decay_card_Dp_KPi0).apply(sel_Dp_KPi0)

# ---------------------------------------------------------------------
# Mode 4: D+ → K_S0 π+,  K_S0 → π+ π-   (≥3 tracks, secondary vertex)
# ---------------------------------------------------------------------
alg_name_Dp_KsPi = "DpToKsPiTag"
alg_Dp_KsPi = Algorithm.new(alg_name_Dp_KsPi)
alg_Dp_KsPi.set_header(["#{alg_name_Dp_KsPi}Alg/#{alg_name_Dp_KsPi}.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                 "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
           .note(:ks0_selection, "K_S0 built from π+π- by a secondary-vertex fit; flight significance " \
                 "L/σ > 2 and |M(π+π-) - m_K_S0| < 12 MeV/c^2 required")

sel_Dp_KsPi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=3"                  # at least three tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:pip, :pim]) {         # K_S0 ← π+π-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :pip]) {    # final fit of the tagged D
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_Dp_KsPi.with_decay_card(decay_card_Dp_KsPi).apply(sel_Dp_KsPi)

# ---------------------------------------------------------------------
# Mode 5: D0 → K_S0 π0,  K_S0 → π+ π-,  π0 → γγ   (≥2 tracks, ≥2 photons)
# ---------------------------------------------------------------------
alg_name_D0_KsPi0 = "D0ToKsPi0Tag"
alg_D0_KsPi0 = Algorithm.new(alg_name_D0_KsPi0)
alg_D0_KsPi0.set_header(["#{alg_name_D0_KsPi0}Alg/#{alg_name_D0_KsPi0}.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })
            .set_alias({ "std::vector<double>" => "Vdouble" })
            .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                  "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
            .note(:ks0_selection, "K_S0 built from π+π- by a secondary-vertex fit; flight significance " \
                  "L/σ > 2 and |M(π+π-) - m_K_S0| < 12 MeV/c^2 required")
            .note(:pi0_mass_window, "reconstructed π0 required in the γγ mass window 0.115-0.150 GeV/c^2")

sel_D0_KsPi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=2"                  # at least two tracks
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:pip, :pim]) {         # K_S0 ← π+π-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # π0 ← γγ, mass constrained
    invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:K_S0, :pi0]) {    # final fit of the tagged D
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_D0_KsPi0.with_decay_card(decay_card_D0_KsPi0).apply(sel_D0_KsPi0)

# ---------------------------------------------------------------------
# Mode 6: D+ → π+ η,  η → γγ   (≥1 track, ≥2 photons)
# ---------------------------------------------------------------------
alg_name_Dp_PiEta = "DpToPiEtaTag"
alg_Dp_PiEta = Algorithm.new(alg_name_Dp_PiEta)
alg_Dp_PiEta.set_header(["#{alg_name_Dp_PiEta}Alg/#{alg_name_Dp_PiEta}.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })
            .set_alias({ "std::vector<double>" => "Vdouble" })
            .note(:background_veto, "cosmic/Bhabha rejection: TOF time difference < 5 ns; " \
                  "e+e- and mu+mu- pairs excluded; require an extra EMC cluster > 50 MeV or an extra MDC track")
            .note(:eta_mass_window, "reconstructed η required in the γγ mass window 0.515-0.575 GeV/c^2")

sel_Dp_PiEta = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"                  # at least one track
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # η ← γγ, mass constrained
    invariant_mass_of(:gamma, :gamma).within(0.515, 0.575)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :eta]) {     # final fit of the tagged D
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_Dp_PiEta.with_decay_card(decay_card_Dp_PiEta).apply(sel_Dp_PiEta)

# =====================================================================
### Execute on datasets (real data + inclusive MC + signal exclusive MC) ###
# mBC and ΔE (±3σ) cuts are applied later at the ROOT level.
# =====================================================================
root_files_D0_KPi   = alg_D0_KPi.execute_on([psi3770_data, psi3770_incMC, exMC_D0_KPi])
root_files_D0_PiPi  = alg_D0_PiPi.execute_on([psi3770_data, psi3770_incMC, exMC_D0_PiPi])
root_files_Dp_KPi0  = alg_Dp_KPi0.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_KPi0])
root_files_Dp_KsPi  = alg_Dp_KsPi.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_KsPi])
root_files_D0_KsPi0 = alg_D0_KsPi0.execute_on([psi3770_data, psi3770_incMC, exMC_D0_KsPi0])
root_files_Dp_PiEta = alg_Dp_PiEta.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_PiEta])