### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # 1.31e9 J/psi events at 3.097 GeV (real data)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # corresponding inclusive MC
cont_3080  = DatasetManager.real_data.find("708_3080")        # 3.08 GeV continuum data (QED-background check)

# --- Decay card for the ST mode: J/psi -> Xi(1530)- anti-Xi+, tag side anti-Xi+ -> anti-Lambda pi+,
#     anti-Lambda -> anti-p pi+; the Xi(1530)- is left undecayed (recoils against the tag) ---
decay_card_st = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi(1530)- anti-Xi+  PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+  PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+  PHSP;
  Enddecay

  End
DECAYCARD

# --- Decay card for the DT radiative mode: Xi(1530)- -> gamma Xi-, Xi- -> Lambda pi-, Lambda -> p pi- ---
decay_card_dt = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi(1530)- anti-Xi+  PHSP;
  Enddecay

  Decay Xi(1530)-
  1.0000 gamma Xi-  PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda0 pi-  PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-  PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+  PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+  PHSP;
  Enddecay

  End
DECAYCARD

# --- 300k-event exclusive MC for each mode ---
exMC_st = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_xi1530_xibar_st"
  config.related_dataset = jpsi_data
  config.events          = 300000
  config.decay_card      = decay_card_st
  config.cross_section   = :default
end

exMC_dt = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_xi1530_xibar_dt_radiative"
  config.related_dataset = jpsi_data
  config.events          = 300000
  config.decay_card      = decay_card_dt
  config.cross_section   = :default
end

### ST selection: single tag anti-Xi+ -> anti-Lambda pi+, anti-Lambda -> anti-p pi+ ###
st_name = "Xi1530ST"
alg_st = Algorithm.new(st_name)
alg_st.set_header(["#{st_name}Alg/#{st_name}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .note(:mass_window, "anti-Lambda(anti-p pi+) mass window +-5 MeV and anti-Xi+(anti-Lambda pi+) mass window +-8 MeV around the nominal masses, applied before the recoil-mass selection")
      .note(:decay_length, "anti-Lambda and anti-Xi+ decay lengths required to be > 0 cm")

st_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"     # at least two pi+
    nChrn ">=1"     # at least one anti-proton
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / anti-p- candidates
    identify :pion,   against: [:kaon, :proton] # pi+ / pi- candidates
    nprm ">=1"      # at least one anti-proton
    npip ">=2"      # at least two pi+
  }
  .remove([:kp <= :chrgp, :km <= :chrgn, :prp <= :chrgp, :prm <= :chrgn])  # remove kaons and protons
  .assign({:chrgp => :pip, :chrgn => :pim})                                 # remaining tracks -> pions
  .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .partial_miss([1]) {                           # Xi(1530)- missing; tag side reconstructed from recoil
    require_recoil_mass 1.44, 1.65               # recoil mass against anti-Xi+ in the Xi(1530) region
  }

alg_st.with_decay_card(decay_card_st).apply(st_selection)
alg_st.execute_on([jpsi_data, jpsi_incMC, exMC_st, cont_3080])

### DT selection (radiative): Xi(1530)- -> gamma Xi-, Xi- -> Lambda pi- ; anti-Xi+ -> anti-Lambda pi+ ###
dt_name = "Xi1530DT"
alg_dt = Algorithm.new(dt_name)
alg_dt.set_header(["#{dt_name}Alg/#{dt_name}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .note(:photon_gap, "photons with 0.80 < |cos(theta)| < 0.86 (EMC barrel-endcap gap) are excluded")
      .note(:decay_length, "Xi- and anti-Xi+ decay lengths required to be > 0 cm")
      .note(:chi2_cut, "4C-fit chi2 < 5 optimised by S/sqrt(S+B) on signal MC vs inclusive MC; only the loose chi2 < 200 is applied in BOSS, the tight cut is applied in the ROOT analysis")

dt_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"     # at least three positive tracks (p, pi+, pi+)
    nChrn ">=3"     # at least three negative tracks (pi-, pi-, anti-p)
  }
  .select_photon {
    tdc_emc_start 0        # EMC timing 0 ...
    tdc_emc_end 14         # ... to 700 ns
    energyThreshold_b 0.025  # barrel deposited energy > 25 MeV
    energyThreshold_e 0.050  # endcap deposited energy > 50 MeV
    angle_to_track 10.0      # angle to nearest charged track > 10 degrees
    nGam ">=1"               # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / anti-p- candidates
    identify :pion,   against: [:kaon, :proton] # pi+ / pi- candidates
    nprp ">=1"
    nprm ">=1"
    npip ">=2"
    npim ">=2"
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])       # remove kaons
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks -> pions
  .secondary_vertex_fit([:prp, :pim]) {          # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :Lambda, :pim, :Lambda_bar, :pip]) {  # 4C fit to gamma + (Lambda pi-) + (anti-Lambda pi+)
    nominal                       # nominal fit, lowest chi2 combination selected automatically
    constrain_four_momentum       # 4C energy-momentum conservation
    invariant_mass_of(:Lambda, :pim).within(1.3137, 1.3297)       # Xi- mass window +-8 MeV
    invariant_mass_of(:Lambda_bar, :pip).within(1.3137, 1.3297)   # anti-Xi+ mass window +-8 MeV
    chi2_cut 200                  # loose in BOSS; tight chi2 < 5 applied in ROOT
  }

alg_dt.with_decay_card(decay_card_dt).apply(dt_selection)
alg_dt.execute_on([jpsi_data, jpsi_incMC, exMC_dt, cont_3080])