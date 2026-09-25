# ============================================================================
#  ψ(3686) → γ χ_c1 ,  χ_c1 → π+π− η_c
#  η_c reconstructed in 16 exclusive decay modes
#  (grouped into charged / K_S / neutral categories)
#  BOSS-side: dataset preparation + event selection (up to the final fit)
# ============================================================================

### ------------------------------ Datasets -------------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data (2712.4×10^6 events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # ψ(3686) inclusive MC

### ------------- η_c channels: 16 modes grouped in 3 categories ------------ ###
# fields: key, category, η_c decay (EvtGen), proton?, kaon?, π0?, η?, fit participants
# every channel shares the χ_c1 side: γ (from ψ) + π+π− (from χ_c1)
channels = [
  ["ppbar",    :charged, "p+ anti-p-",                        true,  false, false, false,
     [:gamma, :pip, :pim, :prp, :prm]],
  ["KpKm2pi",  :charged, "K+ K- pi+ pi-",                     false, true,  false, false,
     [:gamma, :pip, :pim, :kp, :km, :pip, :pim]],
  ["ppbar2pi", :charged, "p+ anti-p- pi+ pi-",                true,  false, false, false,
     [:gamma, :pip, :pim, :prp, :prm, :pip, :pim]],
  ["4pi",      :charged, "pi+ pi- pi+ pi-",                   false, false, false, false,
     [:gamma, :pip, :pim, :pip, :pim, :pip, :pim]],
  ["4K",       :charged, "K+ K- K+ K-",                       false, true,  false, false,
     [:gamma, :pip, :pim, :kp, :km, :kp, :km]],
  ["6pi",      :charged, "pi+ pi- pi+ pi- pi+ pi-",           false, false, false, false,
     [:gamma, :pip, :pim, :pip, :pim, :pip, :pim, :pip, :pim]],
  ["KpKm4pi",  :charged, "K+ K- pi+ pi- pi+ pi-",             false, true,  false, false,
     [:gamma, :pip, :pim, :kp, :km, :pip, :pim, :pip, :pim]],
  ["KSKpi",    :ks,      "K_S0 K+ pi-",                       false, true,  false, false,
     [:gamma, :pip, :pim, :K_S0, :kp, :pim]],
  ["KSK3pi",   :ks,      "K_S0 K+ pi+ pi- pi-",               false, true,  false, false,
     [:gamma, :pip, :pim, :K_S0, :kp, :pip, :pim, :pim]],
  ["KpKmpi0",  :neutral, "K+ K- pi0",                         false, true,  true,  false,
     [:gamma, :pip, :pim, :kp, :km, :pi0]],
  ["ppbarpi0", :neutral, "p+ anti-p- pi0",                    true,  false, true,  false,
     [:gamma, :pip, :pim, :prp, :prm, :pi0]],
  ["KpKmeta",  :neutral, "K+ K- eta",                         false, true,  false, true,
     [:gamma, :pip, :pim, :kp, :km, :eta]],
  ["2pieta",   :neutral, "pi+ pi- eta",                       false, false, false, true,
     [:gamma, :pip, :pim, :pip, :pim, :eta]],
  ["2pi2pi0",  :neutral, "pi+ pi- pi0 pi0",                   false, false, true,  false,
     [:gamma, :pip, :pim, :pip, :pim, :pi0, :pi0]],
  ["4pieta",   :neutral, "pi+ pi- pi+ pi- eta",               false, false, false, true,
     [:gamma, :pip, :pim, :pip, :pim, :pip, :pim, :eta]],
  ["6pipi0",   :neutral, "pi+ pi- pi+ pi- pi+ pi- pi0",       false, false, true,  false,
     [:gamma, :pip, :pim, :pip, :pim, :pip, :pim, :pip, :pim, :pi0]],
]

root_files = []

channels.each do |key, cat, decay, has_proton, has_kaon, needs_pi0, needs_eta, fit_list|

  ### ------------------------ Decay card (EvtGen) ------------------------ ###
  # ψ(3686) → γ χ_c1 ; χ_c1 → π+π− η_c ; η_c → <mode>
  decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 pi+ pi- eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 #{decay} PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
  DECAYCARD

  ### ------------------- Exclusive MC (100k events) ---------------------- ###
  exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gammachic1_etac_#{key}"
    config.related_dataset = psip_data
    config.events          = 100000                     # 100k events per η_c mode
    config.decay_card      = decay_card
    config.cross_section   = :default
  end

  ### ---------------------------- Algorithm ------------------------------ ###
  alg_name = "Chic1Etac_#{key}"
  my_Algorithm = Algorithm.new(alg_name)
  my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })   # √s = 3.686 GeV

  # inexpressible BOSS-side / post-fit procedures
  my_Algorithm
    .note(:background_veto,
          "ψ(3686)→π+π−J/ψ and ψ(3686)→η J/ψ vetoes applied at 3σ; evaluated in ROOT from the stored 4C fit")
    .note(:chic1_selection,
          "χ_c1 selected in the recoil-mass window [3.48, 3.54] GeV/c^2 against the fitted γ")
    .note(:combination_selection,
          "the candidate with the smallest combined χ^2 (4C + 1C + PID + vertex) is retained")
  if cat == :ks
    my_Algorithm.note(:ks_mass_window,
          "K_S0 accepted within |M(π+π−) − m_K_S0| < 12 MeV/c^2 (from the π+π− secondary-vertex fit)")
  end

  # category-dependent nominal χ² cut on the 4C fit
  chi2_cut_value = { charged: 42, ks: 36, neutral: 23 }[cat]

  ### --------------------- Common event selection ------------------------ ###
  event_selection = Selection.new
  event_selection.select_track {          # charged-track selection
        cos_theta 0.93                    # |cosθ| < 0.93
        Vz 10.0                           # |Vz| < 10 cm
        Vr 1.0                            # Vr < 1 cm
        nChrp ">=2"                       # at least 2 positive tracks
        nChrn ">=2"                       # at least 2 negative tracks
      }
    .select_photon {                      # photon selection
        tdc_emc_start 0                   # EMC TDC 0–14
        tdc_emc_end 14
        angle_to_track 10.0               # opening angle to any track > 10°
        energyThreshold_b 0.025           # 25 MeV in the barrel
        energyThreshold_e 0.050           # 50 MeV in the endcap
        nGam ">=1"                        # at least one photon
      }
    .pid(method: :probability) {          # probability-method PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]      # p+ / p̄ vs K and π
        identify :kaon,   against: [:pion, :proton]    # K+ / K− (needed by kaon modes)
        nprp ">=1" if has_proton
        nprm ">=1" if has_proton
      }

  event_selection.remove([:prp <= :chrgp, :prm <= :chrgn])   # drop identified (anti)protons
  event_selection.remove([:kp  <= :chrgp, :km  <= :chrgn])   # drop identified kaons

  # remaining positive / negative tracks → π+ / π−
  event_selection.assign({ :chrgp => :pip, :chrgn => :pim })
                 .remove(:pip) { condition "three_momentum_of(:pip) < 0.4" }  # soft π+ removal
                 .remove(:pim) { condition "three_momentum_of(:pim) < 0.4" }  # soft π− removal

  ### ------------------ K_S0 reconstruction (K_S modes) ------------------ ###
  if cat == :ks
    event_selection.secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
  end

  ### ------------- π0 / η reconstruction (neutral modes, 1C) ------------- ###
  if needs_pi0
    event_selection.kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 200
          npi0 ">=1"
        }
  end
  if needs_eta
    event_selection.kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 200
          neta ">=1"
        }
  end

  ### --------------------- Nominal 4C kinematic fit ---------------------- ###
  event_selection.kinematic_fit(fit_list) {
          nominal                                             # nominal 4C fit
          constrain_four_momentum                             # 4C energy-momentum
          invariant_mass_of(:gamma, :pip, :pim).within(2.80, 3.20)  # γπ+π− window
          chi2_cut chi2_cut_value                             # 42 / 36 / 23 by category
        }
        # competing 4π hypothesis — stores its χ² for the ROOT-level veto
        .kinematic_fit([:pip, :pip, :pim, :pim]) {
          constrain_four_momentum
        }

  ### -------------------- Generate + execute on data --------------------- ###
  my_Algorithm.with_decay_card(decay_card).apply(event_selection)
  root_files << my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])
end