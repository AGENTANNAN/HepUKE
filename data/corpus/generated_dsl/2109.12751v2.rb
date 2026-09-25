# ============================================================================
#  Dataset preparation — 40 BESIII energy points, sqrt(s) = 3.773 - 4.914 GeV
# ============================================================================
data_points = [
  DatasetManager.real_data.find("712_3773"),   # psi(3770)
  DatasetManager.real_data.find("703_4009"),   # 4.009
  DatasetManager.real_data.find("705_4130"),   # 4.130
  DatasetManager.real_data.find("705_4160"),   # 4.160
  DatasetManager.real_data.find("703_4180"),   # 4.178 - 4.600
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4245"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),   # 4.610 - 4.700
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),   # 4.740 - 4.914
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914")
]

# Matching inclusive MC samples (one per energy point)
incMC_points = [
  DatasetManager.inclusive_mc.find("712_3773"),
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4245"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4310"),
  DatasetManager.inclusive_mc.find("705_4315"),
  DatasetManager.inclusive_mc.find("705_4340"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("705_4380"),
  DatasetManager.inclusive_mc.find("703_4390"),
  DatasetManager.inclusive_mc.find("705_4400"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("705_4440"),
  DatasetManager.inclusive_mc.find("703_4470"),
  DatasetManager.inclusive_mc.find("703_4530"),
  DatasetManager.inclusive_mc.find("703_4575"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914")
]

# ============================================================================
#  ConExc decay cards — exclusive continuum (ISR / vacuum polarisation / R)
#  `Particle vpho <ECMS> 0.0` is injected per energy point by the DSL for a
#  multi-point scan and must NOT be declared here.  ConExc is auto-detected
#  from the literal token "ConExc".
# ============================================================================

# Mode 1: e+e- -> K+ K- pi+ pi-   (ConExc built-in mode 14)
decay_card_kkpipi = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 14;
  Enddecay
  End
DECAYCARD

# Mode 2: e+e- -> K+ K- K+ K-     (ConExc built-in mode 16)
decay_card_kkkk = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 16;
  Enddecay
  End
DECAYCARD

# Mode 3: e+e- -> pi+ pi- pi+ pi- (ConExc built-in mode 12)
decay_card_4pi = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 12;
  Enddecay
  End
DECAYCARD

# Mode 4: e+e- -> p pbar pi+ pi-  (form-D, user-supplied cross section)
decay_card_pppipi = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc -1 p+ anti-p- pi+ pi-;
  Enddecay
  End
DECAYCARD

# Mode 5: e+e- -> K+ K- pi+ pi- pi0   (ConExc built-in mode 19)
decay_card_kkpipipi0 = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 19;
  Enddecay
  End
DECAYCARD

# Mode 6: e+e- -> pi+ pi- pi+ pi- pi0 (ConExc built-in mode 17)
decay_card_4pipi0 = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 17;
  Enddecay
  End
DECAYCARD

# Mode 7: e+e- -> p pbar pi+ pi- pi0  (form-D, user-supplied cross section)
decay_card_pppipipi0 = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc -1 p+ anti-p- pi+ pi- pi0;
  Enddecay
  End
DECAYCARD

# Mode 8: e+e- -> K+ K- K+ K- pi0     (form-D, user-supplied cross section)
decay_card_kkkkpi0 = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc -1 K+ K- K+ K- pi0;
  Enddecay
  End
DECAYCARD

# ============================================================================
#  Exclusive ConExc continuum MC: 100k events per mode per energy point
#  (create_exclusive_mc_for returns one ExclusiveMC per energy point)
# ============================================================================
exMC_kkpipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_kkpipi"
  config.events        = 100_000
  config.decay_card    = decay_card_kkpipi
  config.cross_section = :default
end

exMC_kkkk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_kkkk"
  config.events        = 100_000
  config.decay_card    = decay_card_kkkk
  config.cross_section = :default
end

exMC_4pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_4pi"
  config.events        = 100_000
  config.decay_card    = decay_card_4pi
  config.cross_section = :default
end

exMC_pppipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_pppipi"
  config.events        = 100_000
  config.decay_card    = decay_card_pppipi
  config.cross_section = :default
end

exMC_kkpipipi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_kkpipipi0"
  config.events        = 100_000
  config.decay_card    = decay_card_kkpipipi0
  config.cross_section = :default
end

exMC_4pipi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_4pipi0"
  config.events        = 100_000
  config.decay_card    = decay_card_4pipi0
  config.cross_section = :default
end

exMC_pppipipi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_pppipipi0"
  config.events        = 100_000
  config.decay_card    = decay_card_pppipipi0
  config.cross_section = :default
end

exMC_kkkkpi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "conexc_kkkkpi0"
  config.events        = 100_000
  config.decay_card    = decay_card_kkkkpi0
  config.cross_section = :default
end

# ============================================================================
#  Event selection (BOSS) — eight independent final states, one algorithm each
#  Veto windows (GeV/c^2):
#    J/psi  +-15 MeV -> [3.0819, 3.1119]
#    psi2S  +-20 MeV -> [3.6661, 3.7061]
#    D0     +-20 MeV -> [1.8448, 1.8848]
#    chi_c0 +-50 MeV -> [3.3647, 3.4647]
#    K_S0   +-30 MeV -> [0.4676, 0.5276]
# ============================================================================

# ---------------------------------------------------------------------------
# Mode 1: e+e- -> K+ K- pi+ pi-     (4 charged tracks, no photon requirement)
# ---------------------------------------------------------------------------
alg_kkpipi = Algorithm.new("ConExcKKPiPi")
alg_kkpipi.set_header(["ConExcKKPiPiAlg/ConExcKKPiPi.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

sel_kkpipi = Selection.new
  .select_track {
    cos_theta 0.93       # |cos(theta)| < 0.93
    Vz        10.0       # |Vz| < 10 cm
    Vr        1.0        # Vr < 1 cm
    nChrp     "==2"      # exactly 2 positive tracks
    nChrn     "==2"      # exactly 2 negative tracks
    nNet      "==0"      # net charge zero
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+/K-
    identify :pion, against: [:kaon, :proton]   # pi+/pi-
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim]) {      # 4C fit
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :pim).out_of(1.8448, 1.8848)    # |M(K pi) - m_D0|    > 20 MeV
    invariant_mass_of(:kp, :km).out_of(3.3647, 3.4647)     # |M(K K)  - m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)   # |M(pi pi)- m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)   # |M(pi pi)- m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)   # |M(pi pi)- m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_kkpipi.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection on charged tracks) and cos(theta(pi+pi-)) < 0.9")

alg_kkpipi.with_decay_card(decay_card_kkpipi).apply(sel_kkpipi)
alg_kkpipi.execute_on(data_points + incMC_points + exMC_kkpipi)

# ---------------------------------------------------------------------------
# Mode 2: e+e- -> K+ K- K+ K-       (4 charged tracks, no photon requirement)
# ---------------------------------------------------------------------------
alg_kkkk = Algorithm.new("ConExcKKKK")
alg_kkkk.set_header(["ConExcKKKKAlg/ConExcKKKK.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })

sel_kkkk = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+/K-
    nkp  "==2"
    nkm  "==2"
  }
  .kinematic_fit([:kp, :kp, :km, :km]) {       # 4C fit
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).out_of(3.3647, 3.4647)     # |M(K K) - m_chi_c0| > 50 MeV
    chi2_cut 50
  }
alg_kkkk.note(:background_veto,
  "channel-dependent veto E_EMC/p < 0.8 (electron rejection) not expressible in the DSL")

alg_kkkk.with_decay_card(decay_card_kkkk).apply(sel_kkkk)
alg_kkkk.execute_on(data_points + incMC_points + exMC_kkkk)

# ---------------------------------------------------------------------------
# Mode 3: e+e- -> pi+ pi- pi+ pi-   (4 charged tracks, no photon requirement)
# ---------------------------------------------------------------------------
alg_4pi = Algorithm.new("ConExc4Pi")
alg_4pi.set_header(["ConExc4PiAlg/ConExc4Pi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })

sel_4pi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # pi+/pi-
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim]) {   # 4C fit
    nominal
    constrain_four_momentum
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)   # |M(pi pi)- m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pip, :pim, :pim).out_of(3.6661, 3.7061) # |M(4pi) - m_psi(2S)| > 20 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)   # |M(pi pi)- m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)   # |M(pi pi)- m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_4pi.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection) and cos(theta(pi+pi-)) < 0.9")

alg_4pi.with_decay_card(decay_card_4pi).apply(sel_4pi)
alg_4pi.execute_on(data_points + incMC_points + exMC_4pi)

# ---------------------------------------------------------------------------
# Mode 4: e+e- -> p pbar pi+ pi-    (4 charged tracks, no photon requirement)
# ---------------------------------------------------------------------------
alg_pppipi = Algorithm.new("ConExcPPPiPi")
alg_pppipi.set_header(["ConExcPPPiPiAlg/ConExcPPPiPi.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

sel_pppipi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]   # p/pbar
    identify :pion,   against: [:proton]        # pi+/pi-
    nprp "==1"
    nprm "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim]) {   # 4C fit
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :prm).out_of(3.0819, 3.1119)   # |M(p pbar)- m_J/psi| > 15 MeV
    invariant_mass_of(:prp, :prm).out_of(3.6661, 3.7061)   # |M(p pbar)- m_psi(2S)|> 20 MeV
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)   # |M(pi pi) - m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)   # |M(pi pi) - m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)   # |M(pi pi) - m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_pppipi.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection) and cos(theta(pi+pi-)) < 0.9")

alg_pppipi.with_decay_card(decay_card_pppipi).apply(sel_pppipi)
alg_pppipi.execute_on(data_points + incMC_points + exMC_pppipi)

# ---------------------------------------------------------------------------
# Mode 5: e+e- -> K+ K- pi+ pi- pi0   (>=2 photons, 5C fit)
# ---------------------------------------------------------------------------
alg_kkpipipi0 = Algorithm.new("ConExcKKPiPiPi0")
alg_kkpipipi0.set_header(["ConExcKKPiPiPi0Alg/ConExcKKPiPiPi0.h"])
             .set_constant({ "ECMS" => [:double, 3.773] })

sel_kkpipipi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0        # EMC time window 0 - 700 ns
    tdc_emc_end       14
    energyThreshold_b 0.025    # E_gamma > 25 MeV (barrel)
    energyThreshold_e 0.050    # E_gamma > 50 MeV (endcap)
    nGam              ">=2"    # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {   # 5C = 4C + m(pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:kp, :pim).out_of(1.8448, 1.8848)      # |M(K pi) - m_D0|    > 20 MeV
    invariant_mass_of(:kp, :km).out_of(3.3647, 3.4647)       # |M(K K)  - m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)     # |M(pi pi)- m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)     # |M(pi pi)- m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)     # |M(pi pi)- m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_kkpipipi0.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection) and cos(theta(pi+pi-)) < 0.9")

alg_kkpipipi0.with_decay_card(decay_card_kkpipipi0).apply(sel_kkpipipi0)
alg_kkpipipi0.execute_on(data_points + incMC_points + exMC_kkpipipi0)

# ---------------------------------------------------------------------------
# Mode 6: e+e- -> pi+ pi- pi+ pi- pi0   (>=2 photons, 5C fit)
# ---------------------------------------------------------------------------
alg_4pipi0 = Algorithm.new("ConExc4PiPi0")
alg_4pipi0.set_header(["ConExc4PiPi0Alg/ConExc4PiPi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

sel_4pipi0 = Selection.new
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
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim, :gamma, :gamma]) {   # 5C = 4C + m(pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)       # |M(pi pi)- m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pip, :pim, :pim).out_of(3.6661, 3.7061) # |M(4pi) - m_psi(2S)| > 20 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)       # |M(pi pi)- m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)       # |M(pi pi)- m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_4pipi0.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection) and cos(theta(pi+pi-)) < 0.9")

alg_4pipi0.with_decay_card(decay_card_4pipi0).apply(sel_4pipi0)
alg_4pipi0.execute_on(data_points + incMC_points + exMC_4pipi0)

# ---------------------------------------------------------------------------
# Mode 7: e+e- -> p pbar pi+ pi- pi0   (>=2 photons, 5C fit)
# ---------------------------------------------------------------------------
alg_pppipipi0 = Algorithm.new("ConExcPPPiPiPi0")
alg_pppipipi0.set_header(["ConExcPPPiPiPi0Alg/ConExcPPPiPiPi0.h"])
             .set_constant({ "ECMS" => [:double, 3.773] })

sel_pppipipi0 = Selection.new
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
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    identify :pion,   against: [:proton]
    nprp "==1"
    nprm "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma]) {   # 5C = 4C + m(pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:prp, :prm).out_of(3.0819, 3.1119)       # |M(p pbar)- m_J/psi| > 15 MeV
    invariant_mass_of(:prp, :prm).out_of(3.6661, 3.7061)       # |M(p pbar)- m_psi(2S)|> 20 MeV
    invariant_mass_of(:pip, :pim).out_of(3.0819, 3.1119)       # |M(pi pi) - m_J/psi| > 15 MeV
    invariant_mass_of(:pip, :pim).out_of(3.3647, 3.4647)       # |M(pi pi) - m_chi_c0|> 50 MeV
    invariant_mass_of(:pip, :pim).out_of(0.4676, 0.5276)       # |M(pi pi) - m_K_S0|  > 30 MeV
    chi2_cut 50
  }
alg_pppipipi0.note(:background_veto,
  "channel-dependent vetoes not expressible in the DSL: E_EMC/p < 0.8 " \
  "(electron rejection) and cos(theta(pi+pi-)) < 0.9")

alg_pppipipi0.with_decay_card(decay_card_pppipipi0).apply(sel_pppipipi0)
alg_pppipipi0.execute_on(data_points + incMC_points + exMC_pppipipi0)

# ---------------------------------------------------------------------------
# Mode 8: e+e- -> K+ K- K+ K- pi0   (>=2 photons, 5C fit)
# ---------------------------------------------------------------------------
alg_kkkkpi0 = Algorithm.new("ConExcKKKKPi0")
alg_kkkkpi0.set_header(["ConExcKKKKPi0Alg/ConExcKKKKPi0.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })

sel_kkkkpi0 = Selection.new
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
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==2"
    nkm "==2"
  }
  .kinematic_fit([:kp, :kp, :km, :km, :gamma, :gamma]) {   # 5C = 4C + m(pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:kp, :km).out_of(3.3647, 3.4647)     # |M(K K) - m_chi_c0| > 50 MeV
    chi2_cut 50
  }
alg_kkkkpi0.note(:background_veto,
  "channel-dependent veto E_EMC/p < 0.8 (electron rejection) not expressible in the DSL")

alg_kkkkpi0.with_decay_card(decay_card_kkkkpi0).apply(sel_kkkkpi0)
alg_kkkkpi0.execute_on(data_points + incMC_points + exMC_kkkkpi0)