# =====================================================================
# BOSS-side DSL
#   psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi-   / pbar n pi+
#                          and   -> p nbar pi- pi0 / pbar n pi+ pi0
#   for J = 0, 1, 2  ->  12 signal modes; every charge-conjugate channel
#   is given its own selection.
# =====================================================================

### ------------------------------- Datasets ------------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) 3.686 GeV real data (~1.06e8 events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")     # 3.65 GeV continuum real data
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")  # 3.65 GeV continuum inclusive MC

### ---------------------------- Decay cards ---------------------------- ###
# channel A : chi_cJ -> p nbar pi-
card_c0_A = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 P2GC1;
  Enddecay

  Decay chi_c0
  1.0000 p+ anti-n0 pi- PHSP;
  Enddecay

  End
DECAYCARD

card_c1_A = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 p+ anti-n0 pi- PHSP;
  Enddecay

  End
DECAYCARD

card_c2_A = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 P2GC1;
  Enddecay

  Decay chi_c2
  1.0000 p+ anti-n0 pi- PHSP;
  Enddecay

  End
DECAYCARD

# channel B : chi_cJ -> pbar n pi+
card_c0_B = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 P2GC1;
  Enddecay

  Decay chi_c0
  1.0000 anti-p- n0 pi+ PHSP;
  Enddecay

  End
DECAYCARD

card_c1_B = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 anti-p- n0 pi+ PHSP;
  Enddecay

  End
DECAYCARD

card_c2_B = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 P2GC1;
  Enddecay

  Decay chi_c2
  1.0000 anti-p- n0 pi+ PHSP;
  Enddecay

  End
DECAYCARD

# channel C : chi_cJ -> p nbar pi- pi0, pi0 -> gamma gamma
card_c0_C = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 P2GC1;
  Enddecay

  Decay chi_c0
  1.0000 p+ anti-n0 pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

card_c1_C = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 p+ anti-n0 pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

card_c2_C = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 P2GC1;
  Enddecay

  Decay chi_c2
  1.0000 p+ anti-n0 pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# channel D : chi_cJ -> pbar n pi+ pi0, pi0 -> gamma gamma
card_c0_D = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 P2GC1;
  Enddecay

  Decay chi_c0
  1.0000 anti-p- n0 pi+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

card_c1_D = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 anti-p- n0 pi+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

card_c2_D = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 P2GC1;
  Enddecay

  Decay chi_c2
  1.0000 anti-p- n0 pi+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### ------------- Exclusive MC : 500k events for each of the 12 modes ------------- ###
signal_cards = {
  "chi_c0_A" => card_c0_A, "chi_c0_B" => card_c0_B, "chi_c0_C" => card_c0_C, "chi_c0_D" => card_c0_D,
  "chi_c1_A" => card_c1_A, "chi_c1_B" => card_c1_B, "chi_c1_C" => card_c1_C, "chi_c1_D" => card_c1_D,
  "chi_c2_A" => card_c2_A, "chi_c2_B" => card_c2_B, "chi_c2_C" => card_c2_C, "chi_c2_D" => card_c2_D
}

exmc = {}
signal_cards.each do |key, card|
  exmc[key] = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_#{key}"   # 500k-event signal sample per mode
    config.related_dataset = psip_data       # signal is produced on the psi(2S) peak
    config.events          = 500_000
    config.decay_card      = card
    config.cross_section   = :default
  end
end

### ------------------------------ Event selection ------------------------------ ###

# ---- channel A : chi_cJ -> p nbar pi-  (nbar missing, >=1 photon) ----
sel_A = Selection.new
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz        5.0         # |Vz| < 5 cm
    Vr        0.5         # Vr < 0.5 cm
    nChrp     "==1"       # exactly one positive track  (p)
    nChrn     "==1"       # exactly one negative track  (pi-)
    nNet      "==0"       # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0     # at least 10 deg from any charged track
    energyThreshold_b 0.025    # 25 MeV in the barrel
    energyThreshold_e 0.050    # 50 MeV in the endcap
    nGam              ">=1"    # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]    # separate p from K/pi
    identify :pion,   against: [:kaon, :proton]  # separate pi from K/p
    nprp "==1"                                   # exactly one proton candidate
    npim "==1"                                   # exactly one pion candidate
  }
  .select_isolated_photon {
    angle_to_prp_track 30.0   # radiative photon isolated by >30 deg from the proton
    nGam ">=1"
  }
  .for_each(:gamma) do
    where { energy < 0.080 }  # drop photons with E < 80 MeV
    remove
  end
  .kinematic_fit([:gamma, :prp, :n_bar, :pim]) {
    nominal
    miss_track_of :n_bar                  # missing anti-neutron
    constrain_four_momentum               # 1C fit: missing nbar constrained to m_n
    chi2_cut 200                          # loose BOSS-level chi2
  }

# ---- channel B : chi_cJ -> pbar n pi+  (n missing, >=1 photon) ----
sel_B = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        5.0
    Vr        0.5
    nChrp     "==1"       # exactly one positive track  (pi+)
    nChrn     "==1"       # exactly one negative track  (pbar)
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprm "==1"                                   # exactly one anti-proton candidate
    npip "==1"                                   # exactly one pion candidate
  }
  .select_isolated_photon {
    angle_to_prm_track 30.0   # radiative photon isolated by >30 deg from the anti-proton
    nGam ">=1"
  }
  .for_each(:gamma) do
    where { energy < 0.080 }
    remove
  end
  .kinematic_fit([:gamma, :prm, :n, :pip]) {
    nominal
    miss_track_of :n                      # missing neutron
    constrain_four_momentum
    chi2_cut 200
  }

# ---- channel C : chi_cJ -> p nbar pi- pi0  (nbar missing, >=3 photons) ----
sel_C = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        5.0
    Vr        0.5
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"     # 1 radiative photon + 2 photons from pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp "==1"
    npim "==1"
  }
  .select_isolated_photon {
    angle_to_prp_track 30.0
    nGam ">=3"
  }
  .for_each(:gamma) do
    where { energy < 0.080 }
    remove
  end
  # nominal : three-photon hypothesis, all photon combinations tried, smallest chi2 kept
  .kinematic_fit([:gamma, :gamma, :gamma, :prp, :n_bar, :pim]) {
    nominal
    miss_track_of :n_bar
    constrain_four_momentum
    chi2_cut 200
  }
  # competing hypothesis : two photons (no pi0) -> store chi2 for a later ROOT-level veto
  .kinematic_fit([:gamma, :gamma, :prp, :n_bar, :pim]) {
    miss_track_of :n_bar
    constrain_four_momentum
  }
  # competing hypothesis : four photons (extra photon) -> store chi2 for a later veto
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prp, :n_bar, :pim]) {
    miss_track_of :n_bar
    constrain_four_momentum
  }

# ---- channel D : chi_cJ -> pbar n pi+ pi0  (n missing, >=3 photons) ----
sel_D = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        5.0
    Vr        0.5
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprm "==1"
    npip "==1"
  }
  .select_isolated_photon {
    angle_to_prm_track 30.0
    nGam ">=3"
  }
  .for_each(:gamma) do
    where { energy < 0.080 }
    remove
  end
  # nominal : three-photon hypothesis
  .kinematic_fit([:gamma, :gamma, :gamma, :prm, :n, :pip]) {
    nominal
    miss_track_of :n
    constrain_four_momentum
    chi2_cut 200
  }
  # competing hypothesis : two photons
  .kinematic_fit([:gamma, :gamma, :prm, :n, :pip]) {
    miss_track_of :n
    constrain_four_momentum
  }
  # competing hypothesis : four photons
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prm, :n, :pip]) {
    miss_track_of :n
    constrain_four_momentum
  }

### ------------------------------- Algorithms ------------------------------- ###
# chi_c0
alg_c0_A = Algorithm.new("ChiC0PnbarPim")
alg_c0_A.set_header(["ChiC0PnbarPimAlg/ChiC0PnbarPim.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c0_A.with_decay_card(card_c0_A).apply(sel_A.dup)
alg_c0_A.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c0_A"]])

alg_c0_B = Algorithm.new("ChiC0PbarNPip")
alg_c0_B.set_header(["ChiC0PbarNPipAlg/ChiC0PbarNPip.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c0_B.with_decay_card(card_c0_B).apply(sel_B.dup)
alg_c0_B.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c0_B"]])

alg_c0_C = Algorithm.new("ChiC0PnbarPimPi0")
alg_c0_C.set_header(["ChiC0PnbarPimPi0Alg/ChiC0PnbarPimPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c0_C.with_decay_card(card_c0_C).apply(sel_C.dup)
alg_c0_C.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c0_C"]])

alg_c0_D = Algorithm.new("ChiC0PbarNPipPi0")
alg_c0_D.set_header(["ChiC0PbarNPipPi0Alg/ChiC0PbarNPipPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c0_D.with_decay_card(card_c0_D).apply(sel_D.dup)
alg_c0_D.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c0_D"]])

# chi_c1
alg_c1_A = Algorithm.new("ChiC1PnbarPim")
alg_c1_A.set_header(["ChiC1PnbarPimAlg/ChiC1PnbarPim.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c1_A.with_decay_card(card_c1_A).apply(sel_A.dup)
alg_c1_A.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c1_A"]])

alg_c1_B = Algorithm.new("ChiC1PbarNPip")
alg_c1_B.set_header(["ChiC1PbarNPipAlg/ChiC1PbarNPip.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c1_B.with_decay_card(card_c1_B).apply(sel_B.dup)
alg_c1_B.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c1_B"]])

alg_c1_C = Algorithm.new("ChiC1PnbarPimPi0")
alg_c1_C.set_header(["ChiC1PnbarPimPi0Alg/ChiC1PnbarPimPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c1_C.with_decay_card(card_c1_C).apply(sel_C.dup)
alg_c1_C.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c1_C"]])

alg_c1_D = Algorithm.new("ChiC1PbarNPipPi0")
alg_c1_D.set_header(["ChiC1PbarNPipPi0Alg/ChiC1PbarNPipPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c1_D.with_decay_card(card_c1_D).apply(sel_D.dup)
alg_c1_D.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c1_D"]])

# chi_c2
alg_c2_A = Algorithm.new("ChiC2PnbarPim")
alg_c2_A.set_header(["ChiC2PnbarPimAlg/ChiC2PnbarPim.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c2_A.with_decay_card(card_c2_A).apply(sel_A.dup)
alg_c2_A.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c2_A"]])

alg_c2_B = Algorithm.new("ChiC2PbarNPip")
alg_c2_B.set_header(["ChiC2PbarNPipAlg/ChiC2PbarNPip.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c2_B.with_decay_card(card_c2_B).apply(sel_B.dup)
alg_c2_B.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c2_B"]])

alg_c2_C = Algorithm.new("ChiC2PnbarPimPi0")
alg_c2_C.set_header(["ChiC2PnbarPimPi0Alg/ChiC2PnbarPimPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c2_C.with_decay_card(card_c2_C).apply(sel_C.dup)
alg_c2_C.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c2_C"]])

alg_c2_D = Algorithm.new("ChiC2PbarNPipPi0")
alg_c2_D.set_header(["ChiC2PbarNPipPi0Alg/ChiC2PbarNPipPi0.h"]).set_constant({"ECMS" => [:double, 3.686]})
alg_c2_D.with_decay_card(card_c2_D).apply(sel_D.dup)
alg_c2_D.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exmc["chi_c2_D"]])