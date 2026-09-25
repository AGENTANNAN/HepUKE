# BESIII analysis:
#   Ds+ -> f1(1285) e+ nu_e, f1(1285) -> pi+ pi- eta
#   Ds+ -> f1(1420) e+ nu_e, f1(1420) -> K+ K- pi0
# Double-tag technique: reconstruct a Ds- ("single tag") on the tag side, then a
# Ds+ semileptonic decay with a f1 on the signal side plus one transition gamma
# from Ds*+/-.
# Data: 7.33 fb-1 at Ecm in [4.128, 4.226] GeV (8 energy points).

### Datasets ###
data_4130 = DatasetManager.real_data.find("705_4130")   # 4.128 GeV
data_4160 = DatasetManager.real_data.find("705_4160")   # 4.157 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV

incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

all_data = [data_4130, data_4160, data_4180, data_4190,
            data_4200, data_4210, data_4220, data_4230]
all_incMC = [incMC_4130, incMC_4160, incMC_4180, incMC_4190,
             incMC_4200, incMC_4210, incMC_4220, incMC_4230]

### Decay cards ###
# Signal decay card: Ds+ -> f_1(1285) e+ nu_e, f_1(1285) -> pi+ pi- eta
decay_card_f1_1285 = <<~DECAYCARD
  Decay D_s*+
  1.0000  gamma  D_s+                         VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.0000  f_1(1285)  e+  nu_e                 ISGW2;
  Enddecay

  Decay f_1(1285)
  1.0000  pi+  pi-  eta                       PHSP;
  Enddecay

  Decay eta
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Signal decay card: Ds+ -> f_1(1420) e+ nu_e, f_1(1420) -> K+ K- pi0
decay_card_f1_1420 = <<~DECAYCARD
  Decay D_s*+
  1.0000  gamma  D_s+                         VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.0000  f_1(1420)  e+  nu_e                 ISGW2;
  Enddecay

  Decay f_1(1420)
  1.0000  K+  K-  pi0                         PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_f1_1285 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name    = "Dsp_f1_1285_ep_nue"
  config.events         = 500000
  config.decay_card     = decay_card_f1_1285
  config.cross_section  = :default
end

exMC_f1_1420 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name    = "Dsp_f1_1420_ep_nue"
  config.events         = 500000
  config.decay_card     = decay_card_f1_1420
  config.cross_section  = :default
end

exMC_f1_1285.each { |m| m.save_to_config(format: :yaml, file_path: 'exMC_config_f1_1285') }
exMC_f1_1420.each { |m| m.save_to_config(format: :yaml, file_path: 'exMC_config_f1_1420') }

########################################################################
# Algorithm I : Ds+ -> f1(1285) e+ nu_e, f1(1285) -> pi+ pi- eta
########################################################################
alg_f1_1285 = TagAnalysis.new("DsSLf11285")
alg_f1_1285.set_header(["DsSLf11285Alg/DsSLf11285.h"])
           .set_constant({ "ECMS" => [:double, 4.178] })

# Ds- single tag with 12 hadronic modes
alg_f1_1285.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKsK,
          :DstoPiEta,
          :DstoPiEtaPr,
          :DstoKKPiPi0,
          :DstoKsKPiPi,
          :DstoRhoEta,
          :DstoPiEtaPrRhoGam,
          :DstoKsKPi0,
          :DstoPiPiPi,
          :DstoKPiPi,
          :DstoKsKPiPi
  t.charm -1
end

# Signal side: pi+ pi- e+ from Ds+ SL decay, plus eta -> gamma gamma,
# plus one transition photon from Ds*+/- -> gamma Ds+/-
alg_f1_1285.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.photons 3                       # >=3 photons: 2 for eta -> gamma gamma + 1 transition gamma
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e                    # massless neutrino (semileptonic)
end

# 3C kinematic fit: total 4-momentum conservation + eta mass constraint
# on the (gamma gamma) pair from f1(1285) -> pi+ pi- eta.
alg_f1_1285.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
  f.store_fitted_momenta
end

# Inexpressible BOSS-side procedures captured as notes.
alg_f1_1285
  .note(:tag_mrec_window,
        "Per-energy m_rec (system recoiling against tagged Ds-) windows are " \
        "applied to suppress non-Ds*+/- Ds-/+ backgrounds; ranges: " \
        "[2.060,2.150] @4.128, [2.054,2.170] @4.157, [2.048,2.190] @4.178, " \
        "[2.050,2.180] @4.189, [2.046,2.200] @4.199, [2.044,2.210] @4.209, " \
        "[2.042,2.220] @4.219, [2.040,2.220] @4.226 GeV/c^2.")
  .note(:positron_pid,
        "Positron PID uses combined MDC+TOF+EMC probabilities: L'_e > 0 and " \
        "L'_e/(L'_e+L'_pi+L'_K) > 0.8; further require E_EMC/|p_e| > 0.8 c " \
        "to suppress hadron/muon fakes.")
  .note(:pion_momentum_veto,
        "Pion momentum > 100 MeV/c required for non-K_S0 pions to suppress " \
        "D*+/- decay backgrounds.")
  .note(:transition_gamma_selection,
        "The transition photon from Ds*+/- -> gamma Ds+/- is selected via " \
        "a 3C kinematic fit constraining energy-momentum conservation and " \
        "masses of tagged Ds-, signal Ds+, and the parent Ds*+/-. Two " \
        "parent-assignment hypotheses (Ds- as daughter of Ds*- or Ds+ as " \
        "daughter of Ds*+) are tested and the smallest chi2_3C is kept.")
  .note(:mrec2_photon_window,
        "M_rec^2 against photon and tagged Ds- required in (3.82, 4.05) GeV^2/c^4 " \
        "to enforce a real Ds*+/- Ds-/+ event topology.")
  .with_decay_card(decay_card_f1_1285)
  .apply

alg_f1_1285.execute_on(all_data + all_incMC + exMC_f1_1285)

########################################################################
# Algorithm II : Ds+ -> f1(1420) e+ nu_e, f1(1420) -> K+ K- pi0
########################################################################
alg_f1_1420 = TagAnalysis.new("DsSLf11420")
alg_f1_1420.set_header(["DsSLf11420Alg/DsSLf11420.h"])
           .set_constant({ "ECMS" => [:double, 4.178] })

# Ds- single tag with 9 hadronic modes
alg_f1_1420.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,
          :DstoKsK,
          :DstoPiEta,
          :DstoPiEtaPr,
          :DstoKKPiPi0,
          :DstoKsKPiPi,
          :DstoRhoEta,
          :DstoPiEtaPrRhoGam,
          :DstoKsKPi0
  t.charm -1
end

# Signal side: K+ K- e+ from Ds+ SL decay, plus pi0 -> gamma gamma,
# plus one transition photon from Ds*+/- -> gamma Ds+/-
alg_f1_1420.signal_side do |s|
  s.charged(kp: 1, km: 1, ep: 1)
  s.photons 3                       # >=3 photons: 2 for pi0 -> gamma gamma + 1 transition gamma
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e                    # massless neutrino (semileptonic)
end

# 3C kinematic fit: 4-momentum conservation + pi0 mass constraint on gg
alg_f1_1420.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_f1_1420
  .note(:tag_mrec_window,
        "Per-energy m_rec windows against Ds- applied as in f1(1285) mode; " \
        "ranges: [2.060,2.150] @4.128, [2.054,2.170] @4.157, [2.048,2.190] " \
        "@4.178, [2.050,2.180] @4.189, [2.046,2.200] @4.199, [2.044,2.210] " \
        "@4.209, [2.042,2.220] @4.219, [2.040,2.220] @4.226 GeV/c^2.")
  .note(:positron_pid,
        "Positron PID uses combined MDC+TOF+EMC probabilities: L'_e > 0 and " \
        "L'_e/(L'_e+L'_pi+L'_K) > 0.8; further require E_EMC/|p_e| > 0.8 c " \
        "to suppress hadron/muon fakes.")
  .note(:pion_momentum_veto,
        "Pion momentum > 100 MeV/c required for non-K_S0 pions to suppress " \
        "D*+/- decay backgrounds.")
  .note(:transition_gamma_selection,
        "Transition gamma from Ds*+/- -> gamma Ds+/- chosen via 3C kinematic " \
        "fit constraining energy-momentum conservation and Ds-/Ds+/Ds*+/- " \
        "masses; two parent-assignment hypotheses tested, smallest chi2_3C kept.")
  .note(:mrec2_photon_window,
        "M_rec^2 against photon and tagged Ds- required in (3.78, 4.05) GeV^2/c^4.")
  .note(:background_veto,
        "Suppress Ds+ -> phi e+ nu_e (phi -> K+K-) by requiring M(K+K-) > 1.03 " \
        "GeV/c^2. Suppress Ds+ -> K+K-pi+pi0 with pi/e misID by requiring " \
        "M(K+K-pi0 e+) > 1.92 GeV/c^2.")
  .with_decay_card(decay_card_f1_1420)
  .apply

alg_f1_1420.execute_on(all_data + all_incMC + exMC_f1_1420)
