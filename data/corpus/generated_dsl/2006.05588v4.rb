# ==============================================================
# BOSS-side spec — inclusive multiplicity study of ψ(3686) → γ χ_cJ
# (χ_cJ → anything). Multiplicities N_ch / N_sh / N_pi0 are counted
# inclusively; χ_cJ and J/ψ yields are extracted from fits to the
# shower-energy (E_sh) spectrum per multiplicity bin (ROOT-side).
# No 4C kinematic fit is performed.
# ==============================================================

### Datasets ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # inclusive ψ(3686) MC

# Decay card: ψ(2S) → γ χ_c0 with χ_c0 → anything (generic / phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay
    End
DECAYCARD

# 500k-event exclusive MC for the radiative signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamchic0_anything"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Main selection: N_ch >= 1 ###
alg_name = "GamChiCJ"
alg_main = Algorithm.new(alg_name)
alg_main.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_main = Selection.new
  .select_track {                       # charged tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, ≥1 track
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nTot      ">=1"
  }
  .select_photon {                      # neutral showers: EMC TDC [0,14]; 25 MeV barrel / 50 MeV endcap
     tdc_emc_start     0
     tdc_emc_end       14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=1"
  }
  .pid(method: :probability) {          # pions vs kaons/protons, prob > 0.001
     prob_cut 0.001
     identify :pion, against: [:kaon, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # π0 candidates: γγ invariant-mass window
     invariant_mass_of(:gamma, :gamma).between(0.120, 0.145)
  }

alg_main.with_decay_card(decay_card_signal).apply(sel_main)

### N_ch = 0 study: relaxed track requirement + extra background vetoes ###
alg_name0 = "GamChiCJNch0"
alg_nch0 = Algorithm.new(alg_name0)
alg_nch0.set_header(["#{alg_name0}Alg/#{alg_name0}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})
alg_nch0.note(:background_veto,
              "N_ch=0 study: veto events with |Px_neu| > 1.0 GeV/c or |Py_neu| > 1.0 GeV/c (neutral-shower transverse-momentum balance); no direct DSL construct")
       .note(:background_veto_nonpsip,
              "N_ch=0 study: non-ψ(3686) continuum background filter applied on the neutral-multiplicity sample; no direct DSL construct")
       .note(:background_veto_pipi_jpsi,
              "N_ch=0 study: ψ(3686) → π+π- J/ψ background filter applied on the neutral-multiplicity sample; no direct DSL construct")

sel_nch0 = Selection.new
  .select_track {                       # track-quality cuts only; ≥1-track requirement relaxed
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
  }
  .select_photon {
     tdc_emc_start     0
     tdc_emc_end       14
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=1"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :pion, against: [:kaon, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # π0 candidates: γγ invariant-mass window
     invariant_mass_of(:gamma, :gamma).between(0.120, 0.145)
  }

alg_nch0.with_decay_card(decay_card_signal).apply(sel_nch0)

### Execute on the datasets ###
root_files_main = alg_main.execute_on([psip_data, psip_incMC, exMC_signal])
root_files_nch0 = alg_nch0.execute_on([psip_data, psip_incMC, exMC_signal])