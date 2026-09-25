# ============================================================
# 1. Datasets
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(3686) inclusive MC
# continuum (QED) control samples taken off resonance
cont_3080  = DatasetManager.real_data.find("708_3080")      # ~30 pb^-1 @ 3.08 GeV
cont_3650  = DatasetManager.real_data.find("709_3650")      # ~44 pb^-1 @ 3.65 GeV

# ============================================================
# 2. Decay cards (EvtGen format)
# ============================================================
# --- J/psi -> Lambda anti-Lambda ---
decay_card_jpsi_ll = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# --- psi(3686) -> Lambda anti-Lambda ---
decay_card_psip_ll = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# --- J/psi -> Sigma0 anti-Sigma0 ---
decay_card_jpsi_ss = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# --- psi(3686) -> Sigma0 anti-Sigma0 ---
decay_card_psip_ss = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# ============================================================
# 3. Exclusive MC (500k events per channel)
# ============================================================
exMC_jpsi_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ll"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_ll
  config.cross_section   = :default
end

exMC_psip_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ll"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_psip_ll
  config.cross_section   = :default
end

exMC_jpsi_ss = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ss"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi_ss
  config.cross_section   = :default
end

exMC_psip_ss = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ss"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_psip_ss
  config.cross_section   = :default
end

# ============================================================
# 4. Algorithms: one per channel (different final states / beam energies)
# ============================================================
# --- J/psi -> Lambda anti-Lambda ---
alg_name_jll = "JpsiToLLbar"
alg_jpsi_ll = Algorithm.new(alg_name_jll)
alg_jpsi_ll.set_header(["#{alg_name_jll}Alg/#{alg_name_jll}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

# --- psi(3686) -> Lambda anti-Lambda ---
alg_name_pll = "PsipToLLbar"
alg_psip_ll = Algorithm.new(alg_name_pll)
alg_psip_ll.set_header(["#{alg_name_pll}Alg/#{alg_name_pll}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .set_alias({"std::vector<double>" => "Vdouble"})

# --- J/psi -> Sigma0 anti-Sigma0 ---
alg_name_jss = "JpsiToSSbar"
alg_jpsi_ss = Algorithm.new(alg_name_jss)
alg_jpsi_ss.set_header(["#{alg_name_jss}Alg/#{alg_name_jss}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

# --- psi(3686) -> Sigma0 anti-Sigma0 ---
alg_name_pss = "PsipToSSbar"
alg_psip_ss = Algorithm.new(alg_name_pss)
alg_psip_ss.set_header(["#{alg_name_pss}Alg/#{alg_name_pss}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .set_alias({"std::vector<double>" => "Vdouble"})

# ============================================================
# 5. Event selection
# ============================================================
# Common base for the two Lambda anti-Lambda channels (no photon multiplicity required)
ll_selection_base = Selection.new
  .select_track {                 # charged track selection
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz        10.0              # |Vz| < 10 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     ">=2"             # at least two positive tracks
      nChrn     ">=2"             # at least two negative tracks
      nNet      "==0"             # net charge zero
  }
  .select_photon {                # photon selection (no photon multiplicity cut here)
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0      # at least 10 degrees from any charged track
      energyThreshold_b 0.025     # > 25 MeV in the barrel
      energyThreshold_e 0.050     # > 50 MeV in the endcap
  }

# Common base for the two Sigma0 anti-Sigma0 channels (at least two photons)
ss_selection_base = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"
      nChrn     ">=2"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"     # at least two photons for Sigma0 -> gamma Lambda
  }

# ---- Selection chain: psi -> Lambda anti-Lambda ----
jpsi_ll_selection = ll_selection_base.dup
  .pid(method: :probability) {                    # PID: protons against K and pi
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]   # p and anti-p at once
      nprp     ">=1"
      nprm     ">=1"
  }
  .remove([:prp <= :chrgp])                       # remove identified protons
  .remove([:prm <= :chrgn])                       # remove identified anti-protons
  .assign({:chrgp => :pip, :chrgn => :pim})       # remaining tracks treated as pi+/pi-
  .secondary_vertex_fit([:prp, :pim]) {           # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda -> anti-p pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar]) {        # nominal 4C fit on Lambda anti-Lambda
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

psip_ll_selection = jpsi_ll_selection.dup          # identical selection at the psi(3686)

# ---- Selection chain: psi -> Sigma0 anti-Sigma0 ----
jpsi_ss_selection = ss_selection_base.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp     ">=1"
      nprm     ">=1"
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .select_isolated_photon {                       # photons from Sigma0 decays are isolated
      angle_to_prm_track 30.0                     # at least 30 degrees from the anti-proton
      angle_to_prp_track 10.0                     # at least 10 degrees from the proton
      nGam               ">=2"
  }
  .secondary_vertex_fit([:prp, :pim]) {           # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda -> anti-p pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {  # nominal 4C fit on Lambda anti-Lambda gamma gamma
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

psip_ss_selection = jpsi_ss_selection.dup          # identical selection at the psi(3686)

# ============================================================
# 6. Apply selection to the algorithms
# ============================================================
alg_jpsi_ll.with_decay_card(decay_card_jpsi_ll).apply(jpsi_ll_selection)
alg_psip_ll.with_decay_card(decay_card_psip_ll).apply(psip_ll_selection)
alg_jpsi_ss.with_decay_card(decay_card_jpsi_ss).apply(jpsi_ss_selection)
alg_psip_ss.with_decay_card(decay_card_psip_ss).apply(psip_ss_selection)

# ============================================================
# 7. Execute on data + inclusive MC + continuum + exclusive MC
# ============================================================
alg_jpsi_ll.execute_on([jpsi_data, jpsi_incMC, cont_3080, exMC_jpsi_ll])
alg_psip_ll.execute_on([psip_data, psip_incMC, cont_3650, exMC_psip_ll])
alg_jpsi_ss.execute_on([jpsi_data, jpsi_incMC, cont_3080, exMC_jpsi_ss])
alg_psip_ss.execute_on([psip_data, psip_incMC, cont_3650, exMC_psip_ss])