# =============================================================================
# psi(3686) -> pi0 h_c, h_c -> gamma eta_c, eta_c -> 16 exclusive hadronic modes
# BOSS part: dataset preparation + event selection up to the final 4C fit.
# =============================================================================

### ----------------------- Dataset preparation ----------------------- ###
data_3686  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC
data_3650  = DatasetManager.real_data.find("709_3650")     # 3.65 GeV continuum data (~42 pb^-1)
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")  # 3.65 GeV continuum inclusive MC

all_data = [data_3686, incMC_3686, data_3650, incMC_3650]

# Shared top-decay chain (psi(2S) -> pi0 h_c, h_c -> gamma eta_c) plus the
# intermediate-state decays used by the 16 modes (pi0, eta, K_S0).
eta_c_card = lambda do |body|
  <<~CARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    #{body}
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- VSS;
    Enddecay

    End
  CARD
end

# The 16 eta_c decay bodies (final states as listed in the analysis description).
mode_defs = [
  ["pp",      "1.000 p+ anti-p- PHSP;"],
  ["2pipi",   "1.000 pi+ pi- pi+ pi- PHSP;"],
  ["2KK",     "1.000 K+ K- K+ K- PHSP;"],
  ["KKpipi",  "1.000 K+ K- pi+ pi- PHSP;"],
  ["pppipi",  "1.000 p+ anti-p- pi+ pi- PHSP;"],
  ["3pipi",   "1.000 pi+ pi- pi+ pi- pi+ pi- PHSP;"],
  ["KK2pipi", "1.000 K+ K- pi+ pi- pi+ pi- PHSP;"],
  ["KKpi0",   "1.000 K+ K- pi0 PHSP;"],
  ["pppi0",   "1.000 p+ anti-p- pi0 PHSP;"],
  ["KsKpi",   "0.500 K_S0 K+ pi- PHSP;\n0.500 K_S0 K- pi+ PHSP;"],
  ["KsK3pi",  "0.500 K_S0 K+ pi- pi+ pi- PHSP;\n0.500 K_S0 K- pi+ pi+ pi- PHSP;"],
  ["pipiEta", "1.000 pi+ pi- eta PHSP;"],
  ["KKEta",   "1.000 K+ K- eta PHSP;"],
  ["2pipiEta","1.000 pi+ pi- pi+ pi- eta PHSP;"],
  ["pipi2pi0","1.000 pi+ pi- pi0 pi0 PHSP;"],
  ["2pipi2pi0","1.000 pi+ pi- pi+ pi- pi0 pi0 PHSP;"]
]

cards = {}
exMCs = {}
mode_defs.each do |key, body|
  cards[key] = eta_c_card.call(body)
  # One 200k-event exclusive MC sample per eta_c decay mode.
  exMCs[key] = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "etac_#{key}"
    config.related_dataset = data_3686
    config.events          = 200_000
    config.decay_card      = cards[key]
    config.cross_section   = :default
  end
end

### ----------------------- Selection building blocks ----------------------- ###
# Common charged-track + photon selection (mode-dependent multiplicities).
base_sel = lambda do |npos, nneg, ngam|
  Selection.new
    .select_track {
      cos_theta 0.93            # |cos(theta)| < 0.93
      Vz        10.0            # |Vz| < 10 cm
      Vr        1.0             # Vr < 1 cm
      nChrp     "==#{npos}"     # mode-specific positive multiplicity
      nChrn     "==#{nneg}"     # mode-specific negative multiplicity
      nNet      "==0"           # net charge zero
    }
    .select_photon {
      tdc_emc_start     0       # EMC timing window
      tdc_emc_end       14
      angle_to_track    10.0    # >= 10 deg from any charged track
      energyThreshold_b 0.025   # > 25 MeV in the barrel (|cos(theta)| < 0.8)
      energyThreshold_e 0.050   # > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
      nGam              ngam    # mode-dependent photon multiplicity
    }
end

# Common algorithm header/constant setup.
build_alg = lambda do |name|
  a = Algorithm.new(name)
  a.set_header(["#{name}Alg/#{name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})
  a
end

# Attributes handled after the final 4C fit (ROOT-level) -> recorded as notes.
apply_common_notes = lambda do |alg|
  alg.note(:eta_c_mass_window,
           "eta_c invariant-mass window 2.900-3.050 GeV/c^2 is applied in the " \
           "ROOT analysis on the fitted four-momenta; not a BOSS-level cut.")
     .note(:e1_photon_energy,
           "E1 photon (h_c -> gamma eta_c) energy window 0.450-0.550 GeV is " \
           "applied in ROOT on the fitted four-momentum.")
     .note(:chi2_4c_tight_cut,
           "BOSS applies the loose nominal 4C chi2 < 200; the paper's " \
           "channel-dependent chi2_4C (30-70) is applied in the ROOT analysis.")
     .note(:hc_best_candidate,
           "when several h_c candidates survive, the one with the smallest " \
           "combined chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex is " \
           "retained (ROOT-level multi-candidate selection).")
end

### ============================ Mode 1: pp ============================ ###
alg_pp = build_alg.call("EtacPP")
apply_common_notes.call(alg_pp)
sel_pp = base_sel.call(1, 1, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p-
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # soft psi-side pi0 (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :prp, :prm]) {       # 4C: psi(3686) -> pi0 gamma X
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pp.with_decay_card(cards[:pp]).apply(sel_pp)
alg_pp.execute_on(all_data + [exMCs[:pp]])

### ============================ Mode 2: 2(pi+pi-) ============================ ###
alg_2pipi = build_alg.call("Etac2PiPi")
apply_common_notes.call(alg_2pipi)
sel_2pipi = base_sel.call(2, 2, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"
    npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pipi.with_decay_card(cards[:2pipi]).apply(sel_2pipi)
alg_2pipi.execute_on(all_data + [exMCs[:2pipi]])

### ============================ Mode 3: 2(K+K-) ============================ ###
alg_2KK = build_alg.call("Etac2KK")
apply_common_notes.call(alg_2KK)
sel_2KK = base_sel.call(2, 2, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=2"
    nkm ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :kp, :kp, :km, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2KK.with_decay_card(cards[:2KK]).apply(sel_2KK)
alg_2KK.execute_on(all_data + [exMCs[:2KK]])

### ============================ Mode 4: K+K-pi+pi- ============================ ###
alg_KKpipi = build_alg.call("EtacKKPiPi")
apply_common_notes.call(alg_KKpipi)
sel_KKpipi = base_sel.call(2, 2, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KKpipi.with_decay_card(cards[:KKpipi]).apply(sel_KKpipi)
alg_KKpipi.execute_on(all_data + [exMCs[:KKpipi]])

### ============================ Mode 5: pp pi+pi- ============================ ###
alg_pppipi = build_alg.call("EtacPPPiPi")
apply_common_notes.call(alg_pppipi)
sel_pppipi = base_sel.call(2, 2, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
    nprp ">=1"
    nprm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :prp, :prm, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pppipi.with_decay_card(cards[:pppipi]).apply(sel_pppipi)
alg_pppipi.execute_on(all_data + [exMCs[:pppipi]])

### ============================ Mode 6: 3(pi+pi-) ============================ ###
alg_3pipi = build_alg.call("Etac3PiPi")
apply_common_notes.call(alg_3pipi)
sel_3pipi = base_sel.call(3, 3, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"
    npim ">=3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pip, :pip, :pim, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_3pipi.with_decay_card(cards[:3pipi]).apply(sel_3pipi)
alg_3pipi.execute_on(all_data + [exMCs[:3pipi]])

### ============================ Mode 7: K+K-2(pi+pi-) ============================ ###
alg_KK2pipi = build_alg.call("EtacKK2PiPi")
apply_common_notes.call(alg_KK2pipi)
sel_KK2pipi = base_sel.call(3, 3, ">=3")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=2"
    npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :kp, :km, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KK2pipi.with_decay_card(cards[:KK2pipi]).apply(sel_KK2pipi)
alg_KK2pipi.execute_on(all_data + [exMCs[:KK2pipi]])

### ============================ Mode 8: K+K-pi0 ============================ ###
alg_KKpi0 = build_alg.call("EtacKKPi0")
apply_common_notes.call(alg_KKpi0)
sel_KKpi0 = base_sel.call(1, 1, ">=5")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 and eta_c pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :kp, :km, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KKpi0.with_decay_card(cards[:KKpi0]).apply(sel_KKpi0)
alg_KKpi0.execute_on(all_data + [exMCs[:KKpi0]])

### ============================ Mode 9: pp pi0 ============================ ###
alg_pppi0 = build_alg.call("EtacPPPi0")
apply_common_notes.call(alg_pppi0)
sel_pppi0 = base_sel.call(1, 1, ">=5")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :prp, :prm, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pppi0.with_decay_card(cards[:pppi0]).apply(sel_pppi0)
alg_pppi0.execute_on(all_data + [exMCs[:pppi0]])

### ============================ Mode 10: K_S0 K+/- pi-/+ ============================ ###
alg_KsKpi = build_alg.call("EtacKsKPi")
apply_common_notes.call(alg_KsKpi)
alg_KsKpi.note(:ks_secondary_vertex,
               "K_S0 -> pi+ pi- built by a secondary-vertex fit; the K_S0 mass " \
               "window (20 MeV/c^2) and decay-length significance >= 2 sigma are " \
               "applied in the ROOT analysis.")
       .note(:charged_hypothesis_assignment,
             "no PID is applied for this mode; the non-K_S0 positive track is " \
             "assigned to the K+ hypothesis by charge (K_S0 K+ pi-), the " \
             "charge-conjugate K_S0 K- pi+ being handled by the same selection.")
sel_KsKpi = base_sel.call(2, 2, ">=3")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {              # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .assign({:chrgp => :kp})                           # leftover K+ hypothesis
  .kinematic_fit([:pi0, :gamma, :K_S0, :kp, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KsKpi.with_decay_card(cards[:KsKpi]).apply(sel_KsKpi)
alg_KsKpi.execute_on(all_data + [exMCs[:KsKpi]])

### ============================ Mode 11: K_S0 K+/- pi-/+ pi+pi- ============================ ###
alg_KsK3pi = build_alg.call("EtacKsK3Pi")
apply_common_notes.call(alg_KsK3pi)
alg_KsK3pi.note(:ks_secondary_vertex,
                "K_S0 -> pi+ pi- built by a secondary-vertex fit; mass window " \
                "(20 MeV/c^2) and decay length >= 2 sigma applied in ROOT.")
         .note(:charged_hypothesis_assignment,
               "no PID for this mode; among the non-K_S0 charged tracks one K+ " \
               "and one pi+ are assigned by charge/kinematic-fit combination " \
               "(K_S0 K+ pi- pi+ pi- + charge conjugate).")
sel_KsK3pi = base_sel.call(3, 3, ">=3")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {              # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .assign({:chrgp => :kp})                           # leftover K+ hypothesis
  .kinematic_fit([:pi0, :gamma, :K_S0, :kp, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KsK3pi.with_decay_card(cards[:KsK3pi]).apply(sel_KsK3pi)
alg_KsK3pi.execute_on(all_data + [exMCs[:KsK3pi]])

### ============================ Mode 12: pi+pi-eta ============================ ###
alg_pipiEta = build_alg.call("EtacPiPiEta")
apply_common_notes.call(alg_pipiEta)
sel_pipiEta = base_sel.call(1, 1, ">=5")
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # eta -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipiEta.with_decay_card(cards[:pipiEta]).apply(sel_pipiEta)
alg_pipiEta.execute_on(all_data + [exMCs[:pipiEta]])

### ============================ Mode 13: K+K-eta ============================ ###
alg_KKEta = build_alg.call("EtacKKEta")
apply_common_notes.call(alg_KKEta)
sel_KKEta = base_sel.call(1, 1, ">=5")
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # eta -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :kp, :km, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_KKEta.with_decay_card(cards[:KKEta]).apply(sel_KKEta)
alg_KKEta.execute_on(all_data + [exMCs[:KKEta]])

### ============================ Mode 14: 2(pi+pi-)eta ============================ ###
alg_2pipiEta = build_alg.call("Etac2PiPiEta")
apply_common_notes.call(alg_2pipiEta)
sel_2pipiEta = base_sel.call(2, 2, ">=5")
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # eta -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pip, :pim, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pipiEta.with_decay_card(cards[:2pipiEta]).apply(sel_2pipiEta)
alg_2pipiEta.execute_on(all_data + [exMCs[:2pipiEta]])

### ============================ Mode 15: pi+pi-pi0pi0 ============================ ###
alg_pipi2pi0 = build_alg.call("EtacPiPi2Pi0")
apply_common_notes.call(alg_pipi2pi0)
sel_pipi2pi0 = base_sel.call(1, 1, ">=7")
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 + 2 eta_c pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipi2pi0.with_decay_card(cards[:pipi2pi0]).apply(sel_pipi2pi0)
alg_pipi2pi0.execute_on(all_data + [exMCs[:pipi2pi0]])

### ============================ Mode 16: 2(pi+pi-)pi0pi0 ============================ ###
alg_2pipi2pi0 = build_alg.call("Etac2PiPi2Pi0")
apply_common_notes.call(alg_2pipi2pi0)
sel_2pipi2pi0 = base_sel.call(2, 2, ">=7")
  .kalman_kinematic_fit([:gamma, :gamma]) {          # psi-side pi0 + 2 eta_c pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :gamma, :pip, :pip, :pim, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pipi2pi0.with_decay_card(cards[:2pipi2pi0]).apply(sel_2pipi2pi0)
alg_2pipi2pi0.execute_on(all_data + [exMCs[:2pipi2pi0]])