# Study of anti-n p reactions to multi-pion final states
# Novel antineutron source from J/psi -> p pi- anti-n
# Target: proton from beam pipe cooling oil
# Three reactions: anti-n p -> 2pi+ 2pi-, 2pi+ 2pi- pi0, 2pi+ 2pi- 2pi0
# Cross sections measured in 5 momentum intervals (200-1174 MeV/c)
# J/psi data, ECM = 3.097 GeV

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> p pi- anti-n (ST tag for anti-n source)
decay_card_st = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ pi- anti-n PHSP;
    Enddecay

    End
DECAYCARD

# Decay cards for the three signal reactions
# Reaction 1: anti-n p -> pi+ pi+ pi- pi-
decay_card_sig1 = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ pi- anti-n PHSP;
    Enddecay

    End
DECAYCARD

# Reaction 2: anti-n p -> pi+ pi+ pi- pi- pi0
decay_card_sig2 = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ pi- anti-n PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Reaction 3: anti-n p -> pi+ pi+ pi- pi- pi0 pi0
decay_card_sig3 = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ pi- anti-n PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_sig1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_nbarP_to4Pi_3097"
  config.related_dataset = jpsi_data
  config.events = 500_000
  config.decay_card = decay_card_sig1
  config.cross_section = :default
end

exMC_sig2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_nbarP_to5Pi_3097"
  config.related_dataset = jpsi_data
  config.events = 500_000
  config.decay_card = decay_card_sig2
  config.cross_section = :default
end

exMC_sig3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_nbarP_to6Pi_3097"
  config.related_dataset = jpsi_data
  config.events = 500_000
  config.decay_card = decay_card_sig3
  config.cross_section = :default
end

# ---- Algorithm 1: Reaction anti-n p -> 2pi+ 2pi- ----

alg1 = Algorithm.new("nbarPto4Pi")
alg1.set_header(["nbarPto4PiAlg/nbarPto4Pi.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .note(:antineutron_source, "anti-n from J/psi -> p pi- anti-n; anti-n momentum determined from missing mass recoiling against p pi- system")
    .note(:target_proton, "target proton from beam pipe cooling oil; interaction vertex reconstructed from outgoing charged pions")
    .note(:cross_section, "cross section extracted at ROOT level in 5 anti-n momentum intervals from 200 to 1174 MeV/c")

event_selection1 = Selection.new

event_selection1.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=3"
  nChrn ">=3"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nprm ">=0"
  npip ">=2"
  npim ">=2"
}
.remove([:prp <= :chrgp])
.remove([:pip <= :chrgp, :pim <= :chrgn])

# Tag side: reconstruct J/psi -> p pi- + missing anti-n
event_selection1.kinematic_fit([:prp, :pim]) {
  nominal
  constrain_four_momentum
  miss_particle :anti_n
  chi2_cut 200
}

alg1.with_decay_card(decay_card_sig1).apply(event_selection1)
alg1.execute_on([jpsi_data, jpsi_incMC, exMC_sig1])

# ---- Algorithm 2: Reaction anti-n p -> 2pi+ 2pi- pi0 ----

alg2 = Algorithm.new("nbarPto5Pi")
alg2.set_header(["nbarPto5PiAlg/nbarPto5Pi.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .note(:antineutron_source, "anti-n from J/psi -> p pi- anti-n; anti-n momentum determined from missing mass recoiling against p pi- system")
    .note(:target_proton, "target proton from beam pipe cooling oil; interaction vertex reconstructed from outgoing charged pions")
    .note(:cross_section, "cross section extracted at ROOT level in 5 anti-n momentum intervals from 200 to 1174 MeV/c")

event_selection2 = Selection.new

event_selection2.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=3"
  nChrn ">=3"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nprm ">=0"
  npip ">=2"
  npim ">=2"
}
.remove([:prp <= :chrgp])
.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct pi0 -> gamma gamma
event_selection2.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 30
  npi0 ">=1"
}

# Tag side: reconstruct J/psi -> p pi- + missing anti-n
event_selection2.kinematic_fit([:prp, :pim]) {
  nominal
  constrain_four_momentum
  miss_particle :anti_n
  chi2_cut 200
}

alg2.with_decay_card(decay_card_sig2).apply(event_selection2)
alg2.execute_on([jpsi_data, jpsi_incMC, exMC_sig2])

# ---- Algorithm 3: Reaction anti-n p -> 2pi+ 2pi- 2pi0 ----

alg3 = Algorithm.new("nbarPto6Pi")
alg3.set_header(["nbarPto6PiAlg/nbarPto6Pi.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .note(:antineutron_source, "anti-n from J/psi -> p pi- anti-n; anti-n momentum determined from missing mass recoiling against p pi- system")
    .note(:target_proton, "target proton from beam pipe cooling oil; interaction vertex reconstructed from outgoing charged pions")
    .note(:cross_section, "cross section extracted at ROOT level in 5 anti-n momentum intervals from 200 to 1174 MeV/c")

event_selection3 = Selection.new

event_selection3.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=3"
  nChrn ">=3"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nprm ">=0"
  npip ">=2"
  npim ">=2"
}
.remove([:prp <= :chrgp])
.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct pi0 -> gamma gamma (both pi0s)
event_selection3.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 30
  npi0 ">=2"
}

# Tag side: reconstruct J/psi -> p pi- + missing anti-n
event_selection3.kinematic_fit([:prp, :pim]) {
  nominal
  constrain_four_momentum
  miss_particle :anti_n
  chi2_cut 200
}

alg3.with_decay_card(decay_card_sig3).apply(event_selection3)
alg3.execute_on([jpsi_data, jpsi_incMC, exMC_sig3])