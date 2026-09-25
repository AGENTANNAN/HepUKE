# ============================================================
#  Search for baryonic decays of psi(3770) and psi(4040)
#  BOSS part: dataset preparation + event selection
#  (selection chain up to and including the final 4C kinematic fit)
# ============================================================

### ---------------------------- Datasets ---------------------------- ###
psi3770_data = DatasetManager.real_data.find("712_3773")  # psi(3770) data @ 3.773 GeV (2.9 fb^-1)
psi4040_data = DatasetManager.real_data.find("703_4009")  # 4.009 GeV data (psi(4040)), 482 pb^-1
cont_data    = DatasetManager.real_data.find("709_3650")  # continuum data 3.542-3.650 GeV (67 pb^-1), for subtraction

### ------------------------ Decay cards (EvtGen) ------------------------ ###

# ---- psi(3770) -> Lambda Lambda_bar pi+ pi- ----
decay_card_LLpipi_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi+ pi- PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Lambda Lambda_bar pi0 ----
decay_card_LLpi0_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Lambda Lambda_bar eta ----
decay_card_LLeta_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 eta PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Sigma+ anti-Sigma- ----
decay_card_SS_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Sigma0 anti-Sigma0 ----
decay_card_S0S0_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma0 anti-Sigma0 PHSP;
  Enddecay

  Decay Sigma0
  1.0000 Lambda0 gamma PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 anti-Lambda0 gamma PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Xi- anti-Xi+ ----
decay_card_XiXi_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi- anti-Xi+ PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda0 pi- PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# ---- psi(3770) -> Xi0 anti-Xi0 ----
decay_card_Xi0Xi0_3770 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- psi(4040) -> Lambda Lambda_bar pi+ pi-  (extra sample) ----
decay_card_LLpipi_4040 = <<~DECAYCARD
  Decay psi(4040)
  1.0000 Lambda0 anti-Lambda0 pi+ pi- PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# ---- psi(4040) -> Sigma+ anti-Sigma-  (extra sample) ----
decay_card_SS_4040 = <<~DECAYCARD
  Decay psi(4040)
  1.0000 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- psi(4040) -> Xi- anti-Xi+  (extra sample) ----
decay_card_XiXi_4040 = <<~DECAYCARD
  Decay psi(4040)
  1.0000 Xi- anti-Xi+ PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda0 pi- PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

### ------------------------- Exclusive MC samples ------------------------- ###
# 100k events for each of the seven modes at psi(3770)
exMC_LLpipi_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_LLpipi"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_LLpipi_3770
  config.cross_section   = :default
end

exMC_LLpi0_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_LLpi0"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_LLpi0_3770
  config.cross_section   = :default
end

exMC_LLeta_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_LLeta"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_LLeta_3770
  config.cross_section   = :default
end

exMC_SS_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_SS"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_SS_3770
  config.cross_section   = :default
end

exMC_S0S0_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_S0S0"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_S0S0_3770
  config.cross_section   = :default
end

exMC_XiXi_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_XiXi"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_XiXi_3770
  config.cross_section   = :default
end

exMC_Xi0Xi0_3770 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi3770_Xi0Xi0"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_Xi0Xi0_3770
  config.cross_section   = :default
end

# Extra 100k-event psi(4040) samples for Lambda Lambda_bar pi+ pi-, Sigma+ anti-Sigma- and Xi- anti-Xi+
exMC_LLpipi_4040 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi4040_LLpipi"
  config.related_dataset = psi4040_data
  config.events          = 100_000
  config.decay_card      = decay_card_LLpipi_4040
  config.cross_section   = :default
end

exMC_SS_4040 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi4040_SS"
  config.related_dataset = psi4040_data
  config.events          = 100_000
  config.decay_card      = decay_card_SS_4040
  config.cross_section   = :default
end

exMC_XiXi_4040 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_psi4040_XiXi"
  config.related_dataset = psi4040_data
  config.events          = 100_000
  config.decay_card      = decay_card_XiXi_4040
  config.cross_section   = :default
end

### ============================ Event selection ============================ ###

# ------------------------------------------------------------------
# Mode 1 : psi(3770)/psi(4040) -> Lambda Lambda_bar pi+ pi-
# ------------------------------------------------------------------
alg_LLpipi = Algorithm.new("LLpipi")
alg_LLpipi.set_header(["LLpipiAlg/LLpipi.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:ecms_per_dataset,
                "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
                "continuum (3.542-3.650 GeV) data sets are processed by the same selection " \
                "and require their own beam energy in the 4C fit")

sel_LLpipi = Selection.new
  .select_track {
    cos_theta 0.93   # |cos(theta)| < 0.93
    Vz        10.0   # |Vz| < 10 cm
    Vr        1.0    # Vr < 1 cm
    nChrp     "==3"  # (3+,3-) for Lambda Lambda_bar pi+ pi-
    nChrn     "==3"
    nNet      "==0"  # net charge zero
  }
  .pid(method: :probability) {
    prob_cut 0.001                                        # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]             # p+ / anti-p- vs K and pi
    nprp ">=1"                                            # at least one proton
    nprm ">=1"                                            # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])               # remove identified (anti-)protons
  .assign({:chrgp => :pip, :chrgn => :pim})               # remaining tracks -> pi+ / pi-
  .secondary_vertex_fit([:prp, :pim]) {                   # Lambda -> p pi- (mass-difference minimization)
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {                   # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {    # 4C fit [Lambda Lambda_bar pi+ pi-]
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_LLpipi.with_decay_card(decay_card_LLpipi_3770).apply(sel_LLpipi)

# ------------------------------------------------------------------
# Mode 2 : psi(3770) -> Lambda Lambda_bar pi0
# ------------------------------------------------------------------
alg_LLpi0 = Algorithm.new("LLpi0")
alg_LLpi0.set_header(["LLpi0Alg/LLpi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_per_dataset,
               "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
               "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_LLpi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"   # (2+,2-)
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0   # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025  # barrel  E > 25 MeV
    energyThreshold_e 0.050  # endcap  E > 50 MeV
    nGam              ">=2"  # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {              # 1C pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) {         # 4C fit [Lambda Lambda_bar pi0]
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_LLpi0.with_decay_card(decay_card_LLpi0_3770).apply(sel_LLpi0)

# ------------------------------------------------------------------
# Mode 3 : psi(3770) -> Lambda Lambda_bar eta
# ------------------------------------------------------------------
alg_LLeta = Algorithm.new("LLeta")
alg_LLeta.set_header(["LLetaAlg/LLeta.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:ecms_per_dataset,
               "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
               "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_LLeta = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
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
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {              # 1C eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :eta]) {         # 4C fit [Lambda Lambda_bar eta]
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_LLeta.with_decay_card(decay_card_LLeta_3770).apply(sel_LLeta)

# ------------------------------------------------------------------
# Mode 4 : psi(3770) -> Sigma+ anti-Sigma-  (final state p anti-p pi0 pi0)
# ------------------------------------------------------------------
alg_SS = Algorithm.new("SigmaSigma")
alg_SS.set_header(["SigmaSigmaAlg/SigmaSigma.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:ecms_per_dataset,
            "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
            "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_SS = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"   # (1+,1-)
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"  # Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0  => >= 4 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"        # exactly one proton
    nprm "==1"        # exactly one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .kalman_kinematic_fit([:gamma, :gamma]) {              # 1C pi0 -> gamma gamma (two pi0 candidates)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # No Lambda / anti-Lambda secondary vertex fit in this mode
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {             # 4C fit [p anti-p pi0 pi0]
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_SS.with_decay_card(decay_card_SS_3770).apply(sel_SS)

# ------------------------------------------------------------------
# Mode 5 : psi(3770) -> Sigma0 anti-Sigma0  (final state Lambda Lambda_bar gamma gamma)
# ------------------------------------------------------------------
alg_S0S0 = Algorithm.new("Sigma0Sigma0")
alg_S0S0.set_header(["Sigma0Sigma0Alg/Sigma0Sigma0.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:ecms_per_dataset,
              "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
              "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_S0S0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"   # (2+,2-)
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"  # Sigma0 -> Lambda gamma (x2)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {  # 4C fit [Lambda Lambda_bar gamma gamma]
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:Lambda).within(1.111, 1.121)         # Lambda mass window
    invariant_mass_of(:Lambda_bar).within(1.111, 1.121)     # anti-Lambda mass window
  }

alg_S0S0.with_decay_card(decay_card_S0S0_3770).apply(sel_S0S0)

# ------------------------------------------------------------------
# Mode 6 : psi(3770) -> Xi- anti-Xi+  (final state Lambda Lambda_bar pi+ pi- pi- pi+)
# ------------------------------------------------------------------
alg_XiXi = Algorithm.new("XiXi")
alg_XiXi.set_header(["XiXiAlg/XiXi.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:ecms_per_dataset,
              "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
              "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_XiXi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==3"   # (3+,3-)
    nChrn     "==3"
    nNet      "==0"
  }
  # no photon selection: no photons in this final state
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {       # 4C fit [Lambda Lambda_bar pi+ pi-]
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:Lambda).within(1.111, 1.121)          # Lambda mass window
    invariant_mass_of(:Lambda_bar).within(1.111, 1.121)      # anti-Lambda mass window
    # the fit automatically keeps the smallest-chi2 combination (best Xi- / anti-Xi+ assignment)
  }

alg_XiXi.with_decay_card(decay_card_XiXi_3770).apply(sel_XiXi)

# ------------------------------------------------------------------
# Mode 7 : psi(3770) -> Xi0 anti-Xi0  (final state Lambda Lambda_bar pi0 pi0)
# ------------------------------------------------------------------
alg_Xi0Xi0 = Algorithm.new("Xi0Xi0")
alg_Xi0Xi0.set_header(["Xi0Xi0Alg/Xi0Xi0.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:ecms_per_dataset,
                "ECMS carries the psi(3770) c.m. energy; the psi(4040) (4.009 GeV) and " \
                "continuum (3.542-3.650 GeV) data sets require their own beam energy in the 4C fit")

sel_Xi0Xi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"   # (2+,2-)
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"  # Xi0 -> Lambda pi0 (x2)  => >= 4 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {                  # 1C pi0 -> gamma gamma (two pi0 candidates)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {       # 4C fit [Lambda Lambda_bar pi0 pi0]
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:Lambda).within(1.111, 1.121)          # Lambda mass window
    invariant_mass_of(:Lambda_bar).within(1.111, 1.121)      # anti-Lambda mass window
  }

alg_Xi0Xi0.with_decay_card(decay_card_Xi0Xi0_3770).apply(sel_Xi0Xi0)

### ---------------- Run each selection on data (+ MC) ---------------- ###
# Each selection is applied to the psi(3770) data, the 4.009 GeV data and the
# continuum data; the continuum sample is used for background subtraction.

root_files_LLpipi = alg_LLpipi.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_LLpipi_3770, exMC_LLpipi_4040])

root_files_LLpi0 = alg_LLpi0.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_LLpi0_3770])

root_files_LLeta = alg_LLeta.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_LLeta_3770])

root_files_SS = alg_SS.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_SS_3770, exMC_SS_4040])

root_files_S0S0 = alg_S0S0.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_S0S0_3770])

root_files_XiXi = alg_XiXi.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_XiXi_3770, exMC_XiXi_4040])

root_files_Xi0Xi0 = alg_Xi0Xi0.execute_on(
  [psi3770_data, psi4040_data, cont_data, exMC_Xi0Xi0_3770])