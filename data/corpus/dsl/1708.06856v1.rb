# Search for D+ -> D0 e+ nu_e at BESIII (2.93 fb-1 at sqrt(s)=3.773 GeV)
# Double-tag technique: D- tag on one side, D0 reconstructed on signal side;
# positron and neutrino are not reconstructed (positron momentum < 5 MeV/c).

### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for D+ -> D0 e+ nu_e with D0 -> K- pi+
decay_card_Kpi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e  PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+  PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi-        PHSP;
    0.1667 K+ pi- pi- pi0    PHSP;
    0.1667 K_S0 pi-          PHSP;
    0.1667 K_S0 pi- pi0      PHSP;
    0.1667 K_S0 pi+ pi- pi-  PHSP;
    0.1667 K+ K- pi-         PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D+ -> D0 e+ nu_e with D0 -> K- pi+ pi+ pi-
decay_card_K3pi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e  PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi+ pi-  PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi-        PHSP;
    0.1667 K+ pi- pi- pi0    PHSP;
    0.1667 K_S0 pi-          PHSP;
    0.1667 K_S0 pi- pi0      PHSP;
    0.1667 K_S0 pi+ pi- pi-  PHSP;
    0.1667 K+ K- pi-         PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D+ -> D0 e+ nu_e with D0 -> K- pi+ pi0
decay_card_Kpipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-  PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e  PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi0  PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi-        PHSP;
    0.1667 K+ pi- pi- pi0    PHSP;
    0.1667 K_S0 pi-          PHSP;
    0.1667 K_S0 pi- pi0      PHSP;
    0.1667 K_S0 pi+ pi- pi-  PHSP;
    0.1667 K+ K- pi-         PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_Kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DpToD0enue_D0toKpi"
  config.related_dataset = psi3770_data
  config.events          = 200000
  config.decay_card      = decay_card_Kpi
  config.cross_section   = :default
end

exMC_K3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DpToD0enue_D0toK3pi"
  config.related_dataset = psi3770_data
  config.events          = 200000
  config.decay_card      = decay_card_K3pi
  config.cross_section   = :default
end

exMC_Kpipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DpToD0enue_D0toKpipi0"
  config.related_dataset = psi3770_data
  config.events          = 200000
  config.decay_card      = decay_card_Kpipi0
  config.cross_section   = :default
end

### Tag-based event selection ###

# Six D- tag modes (charm = -1): K+pi-pi-, K+pi-pi-pi0, KS pi-, KS pi- pi0,
# KS pi+ pi- pi-, K+ K- pi-. All modes are DTagAlg standard channels.
d_tag_modes = [
  :DptoKPiPi,        # K- pi+ pi+ (charge conjugate: K+ pi- pi-)
  :DptoKPiPiPi0,     # K- pi+ pi+ pi0
  :DptoKsPi,         # K_S0 pi+
  :DptoKsPiPi0,      # K_S0 pi+ pi0
  :DptoKsPiPiPi,     # K_S0 pi+ pi+ pi-
  :DptoKKPi          # K+ K- pi+
]

# ---- Mode 1: D0 -> K- pi+ ----
alg_Kpi = TagAnalysis.new("DpToD0enueKpi")
alg_Kpi.set_header(["DpToD0enueKpiAlg/DpToD0enueKpi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })

alg_Kpi.tag_side(:Dplus) do |t|
  t.modes(*d_tag_modes)
  t.charm -1
end

alg_Kpi.signal_side do |s|
  s.charged(km: 1, pip: 1)      # D0 -> K- pi+
  s.require_charge 0
  s.missing :nu_e               # e+ (soft, <5 MeV/c) and nu_e represented as one massless missing
  s.min_photon_angle 10.0
end

alg_Kpi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_Kpi.note(:signal_optimization,
  "Optimized probability requirement > 0.37 on the product of normalized D0 momentum and " \
  "normalized observed D- D0 energy distributions; applied on ROOT side.")
   .note(:tag_mbc_selection,
  "ST D- candidates require MBC > 1.83 GeV/c^2 and mode-dependent ~3 sigma DeltaE windows; " \
  "smallest |DeltaE| chosen per mode. Applied on ROOT side (store-not-cut in BOSS).")
   .note(:positron_not_reconstructed,
  "Positron from D+ -> D0 e+ nu_e has momentum < 5 MeV/c and is not reconstructed; " \
  "treated together with nu_e as a single massless missing four-vector.")
   .with_decay_card(decay_card_Kpi)
   .apply

alg_Kpi.execute_on([psi3770_data, psi3770_incMC, exMC_Kpi])


# ---- Mode 2: D0 -> K- pi+ pi+ pi- ----
alg_K3pi = TagAnalysis.new("DpToD0enueK3pi")
alg_K3pi.set_header(["DpToD0enueK3piAlg/DpToD0enueK3pi.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })

alg_K3pi.tag_side(:Dplus) do |t|
  t.modes(*d_tag_modes)
  t.charm -1
end

alg_K3pi.signal_side do |s|
  s.charged(km: 1, pip: 2, pim: 1)   # D0 -> K- pi+ pi+ pi-
  s.require_charge 0
  s.missing :nu_e
  s.min_photon_angle 10.0
end

alg_K3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_K3pi.note(:signal_optimization,
  "Optimized probability requirement > 0.34 on the product of normalized D0 momentum and " \
  "normalized observed D- D0 energy distributions; applied on ROOT side.")
   .note(:tag_mbc_selection,
  "ST D- candidates require MBC > 1.83 GeV/c^2 and mode-dependent ~3 sigma DeltaE windows; " \
  "smallest |DeltaE| chosen per mode. Applied on ROOT side.")
   .note(:positron_not_reconstructed,
  "Positron and nu_e combined into one massless missing four-vector (see Mode 1 note).")
   .with_decay_card(decay_card_K3pi)
   .apply

alg_K3pi.execute_on([psi3770_data, psi3770_incMC, exMC_K3pi])


# ---- Mode 3: D0 -> K- pi+ pi0 ----
alg_Kpipi0 = TagAnalysis.new("DpToD0enueKpipi0")
alg_Kpipi0.set_header(["DpToD0enueKpipi0Alg/DpToD0enueKpipi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

alg_Kpipi0.tag_side(:Dplus) do |t|
  t.modes(*d_tag_modes)
  t.charm -1
end

alg_Kpipi0.signal_side do |s|
  s.photons 2                      # two photons for pi0 -> gamma gamma
  s.charged(km: 1, pip: 1)         # D0 -> K- pi+ (pi0)
  s.require_charge 0
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_Kpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:km, :pip, :gamma, :gamma).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_Kpipi0.note(:signal_optimization,
  "Optimized probability requirement > 0.54 on the product of normalized D0 momentum and " \
  "normalized observed D- D0 energy distributions; applied on ROOT side.")
   .note(:tag_mbc_selection,
  "ST D- candidates require MBC > 1.83 GeV/c^2 and mode-dependent ~3 sigma DeltaE windows; " \
  "smallest |DeltaE| chosen per mode. Applied on ROOT side.")
   .note(:pi0_selection,
  "pi0 candidates: diphoton mass in (0.110, 0.155) GeV/c^2; both-photons-in-endcap rejected; " \
  "1C kinematic fit with chi^2 < 20.")
   .note(:positron_not_reconstructed,
  "Positron and nu_e combined into one massless missing four-vector (see Mode 1 note).")
   .with_decay_card(decay_card_Kpipi0)
   .apply

alg_Kpipi0.execute_on([psi3770_data, psi3770_incMC, exMC_Kpipi0])
