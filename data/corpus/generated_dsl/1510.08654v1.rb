# frozen_string_literal: true
# e+e- -> gamma_ISR/FSR mu+mu- : measurement of E_cms from 3.810 to 4.600 GeV
# BOSS (dataset preparation + event selection) part only.

### ===================== Dataset description ===================== ###
# 18 nominal real-data samples spanning 3.810 - 4.600 GeV (25 physical data sets;
# the separately acquired sub-samples at 4.009, 4.230, 4.260 and 4.420 GeV are
# merged automatically by using the sample name).
data_points = %w[
  703_3810 703_3900 703_4009 703_4090 703_4180 703_4190
  703_4200 703_4210 703_4220 703_4230 703_4237 703_4246
  703_4260 703_4270 703_4280 703_4360 703_4420 703_4600
].map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC for the 4009, 4190scan, 4210scan, 4220scan, 4230, 4260, 4360, 4420 and 4600 points
incMC_points = %w[
  703_4009 703_4190scan 703_4210scan 703_4220scan 703_4230
  703_4260 703_4360 703_4420 703_4600
].map { |name| DatasetManager.inclusive_mc.find(name) }

### ===================== Decay cards ===================== ###
# Di-muon signal, ISR/FSR switched ON (radiative corrections to the VLL decay)
decay_card_dimuon_radon = <<~DECAYCARD
  Decay psi(4260)
  1.0 mu+ mu- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Di-muon signal, ISR/FSR switched OFF
decay_card_dimuon_radoff = <<~DECAYCARD
  Decay psi(4260)
  1.0 mu+ mu- VLL;
  Enddecay
  End
DECAYCARD

# J/psi momentum check: e+e- -> gamma_ISR J/psi, J/psi -> mu+mu-  (FSR ON)
decay_card_jpsi_fsron = <<~DECAYCARD
  Decay psi(4260)
  1.0 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0 mu+ mu- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# J/psi momentum check: e+e- -> gamma_ISR J/psi, J/psi -> mu+mu-  (FSR OFF)
decay_card_jpsi_fsroff = <<~DECAYCARD
  Decay psi(4260)
  1.0 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0 mu+ mu- VLL;
  Enddecay
  End
DECAYCARD

# D0 -> K- pi+ / anti-D0 -> K+ pi-
decay_card_d0kpi = <<~DECAYCARD
  Decay psi(4260)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 K- pi+ PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# e+e- -> pi+ pi- K+ K-
decay_card_pipikk = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- K+ K- PHSP;
  Enddecay
  End
DECAYCARD

# e+e- -> pi+ pi- p anti-p
decay_card_pipippbar = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- p+ anti-p- PHSP;
  Enddecay
  End
DECAYCARD

### ===================== Exclusive MC (50k events per energy point, per mode) ===================== ###
exMC_dimuon_radon = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_dimuon_isrfsr_on"
  config.events        = 50_000
  config.decay_card    = decay_card_dimuon_radon
  config.cross_section = :default
end

exMC_dimuon_radoff = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_dimuon_isrfsr_off"
  config.events        = 50_000
  config.decay_card    = decay_card_dimuon_radoff
  config.cross_section = :default
end

exMC_jpsi_fsron = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_jpsi_fsr_on"
  config.events        = 50_000
  config.decay_card    = decay_card_jpsi_fsron
  config.cross_section = :default
end

exMC_jpsi_fsroff = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_jpsi_fsr_off"
  config.events        = 50_000
  config.decay_card    = decay_card_jpsi_fsroff
  config.cross_section = :default
end

exMC_d0kpi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_d0_to_kpi"
  config.events        = 50_000
  config.decay_card    = decay_card_d0kpi
  config.cross_section = :default
end

exMC_pipikk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pipi_kk"
  config.events        = 50_000
  config.decay_card    = decay_card_pipikk
  config.cross_section = :default
end

exMC_pipippbar = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pipi_ppbar"
  config.events        = 50_000
  config.decay_card    = decay_card_pipippbar
  config.cross_section = :default
end

### ===================== Event selection (BOSS) ===================== ###

# ------------------------------------------------
# Channel 1 : e+e- -> gamma_ISR/FSR mu+mu-
# ------------------------------------------------
alg_dimuon = Algorithm.new("ECmsDimuon")
alg_dimuon.set_header(["ECmsDimuonAlg/ECmsDimuon.h"])
          # NOTE: the scan runs on 18 points; ECMS is fixed here as a single constant,
          # the per-point beam energy is resolved from the dataset at run time.
          .set_constant({ "ECMS" => [:double, 3.810] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .note(:energy_scan, "Measurement of E_cms from e+e- -> gamma_ISR/FSR mu+mu- over 18 data points
            (3.810 - 4.600 GeV); the algorithm package is executed once per energy point and the
            c.m. energy is taken from the dataset, not from the single ECMS constant.")

sel_dimuon = Selection.new
sel_dimuon.select_track {
            cos_theta 0.80      # |cos(theta)| < 0.80
            Vz        10.0      # |Vz| < 10 cm
            Vr        10.0      # Vr < 1 cm (10 mm)
            nChrp     "==1"     # exactly one positive track
            nChrn     "==1"     # exactly one negative track
            nNet      "==0"     # net charge zero
          }
          # per-track EMC energy deposition < 0.4 GeV (remove tracks with a hard shower)
          .for_each(:chrgp) {
            where { Eraw > 0.4 }
            remove
          }
          .for_each(:chrgn) {
            where { Eraw > 0.4 }
            remove
          }
          .pid(method: :probability) {
            # high-momentum lepton pair: p > 1.0 -> lepton; EMC energy > 0.6 -> electron, else muon
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            nlp "==1"           # exactly one l+
            nlm "==1"           # exactly one l-
          }
          # nominal 4C fit : final-state four-momenta constrained to the c.m. four-momentum
          .kinematic_fit([:lp, :lm]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_dimuon.with_decay_card(decay_card_dimuon_radon).apply(sel_dimuon)
alg_dimuon.execute_on(data_points + incMC_points + exMC_dimuon_radon + exMC_dimuon_radoff)

# ------------------------------------------------
# Channel 2 : e+e- -> gamma_ISR J/psi, J/psi -> mu+mu-   (J/psi momentum check)
# ------------------------------------------------
alg_jpsi = Algorithm.new("ECmsJPsi")
alg_jpsi.set_header(["ECmsJPsiAlg/ECmsJPsi.h"])
        .set_constant({ "ECMS" => [:double, 3.810] })
        .set_alias({ "std::vector<double>" => "Vdouble" })

sel_jpsi = Selection.new
sel_jpsi.select_track {
           cos_theta 0.80       # |cos(theta)| < 0.80
           Vz        10.0       # |Vz| < 10 cm
           Vr        10.0       # Vr < 1 cm (10 mm)
           nChrp     "==1"      # exactly one positive track
           nChrn     "==1"      # exactly one negative track
           nNet      "==0"      # net charge zero
         }
         .pid(method: :probability) {
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .kinematic_fit([:lp, :lm]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }

alg_jpsi.with_decay_card(decay_card_jpsi_fsron).apply(sel_jpsi)
alg_jpsi.execute_on(data_points + incMC_points + exMC_jpsi_fsron + exMC_jpsi_fsroff)

# ------------------------------------------------
# Channel 3 : D0 -> K- pi+ / anti-D0 -> K+ pi-  (validation)
# ------------------------------------------------
alg_d0kpi = Algorithm.new("D0KPiCheck")
alg_d0kpi.set_header(["D0KPiCheckAlg/D0KPiCheck.h"])
         .set_constant({ "ECMS" => [:double, 3.810] })
         .set_alias({ "std::vector<double>" => "Vdouble" })

sel_d0kpi = Selection.new
sel_d0kpi.select_track {
           cos_theta 0.93       # |cos(theta)| < 0.93
           Vz        10.0       # |Vz| < 10 cm
           Vr        10.0       # Vr < 1 cm (10 mm)
           nChrp     "==1"      # one positive track
           nChrn     "==1"      # one negative track
           nNet      "==0"      # net charge zero
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]   # K+ and K-
           identify :pion, against: [:kaon, :proton]   # pi+ and pi-
           nkm  "==1"           # exactly one K-
           npip "==1"           # exactly one pi+   (D0 -> K- pi+)
         }
         .kinematic_fit([:km, :pip]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }

alg_d0kpi.with_decay_card(decay_card_d0kpi).apply(sel_d0kpi)
alg_d0kpi.execute_on(data_points + incMC_points + exMC_d0kpi)

# ------------------------------------------------
# Channel 4 : e+e- -> pi+ pi- K+ K-  (validation)
# ------------------------------------------------
alg_pipikk = Algorithm.new("PiPiKKCheck")
alg_pipikk.set_header(["PiPiKKCheckAlg/PiPiKKCheck.h"])
          .set_constant({ "ECMS" => [:double, 3.810] })
          .set_alias({ "std::vector<double>" => "Vdouble" })

sel_pipikk = Selection.new
sel_pipikk.select_track {
            cos_theta 0.93      # |cos(theta)| < 0.93
            Vz        10.0      # |Vz| < 10 cm
            Vr        10.0      # Vr < 1 cm (10 mm)
            nChrp     "==2"     # two positive tracks
            nChrn     "==2"     # two negative tracks
            nNet      "==0"     # net charge zero
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
            nkp  "==1"          # exactly one K+
            nkm  "==1"          # exactly one K-
            npip "==1"          # exactly one pi+
            npim "==1"          # exactly one pi-
          }
          .kinematic_fit([:kp, :km, :pip, :pim]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_pipikk.with_decay_card(decay_card_pipikk).apply(sel_pipikk)
alg_pipikk.execute_on(data_points + incMC_points + exMC_pipikk)

# ------------------------------------------------
# Channel 5 : e+e- -> pi+ pi- p anti-p  (validation)
# ------------------------------------------------
alg_pipippbar = Algorithm.new("PiPiPPbarCheck")
alg_pipippbar.set_header(["PiPiPPbarCheckAlg/PiPiPPbarCheck.h"])
             .set_constant({ "ECMS" => [:double, 3.810] })
             .set_alias({ "std::vector<double>" => "Vdouble" })

sel_pipippbar = Selection.new
sel_pipippbar.select_track {
               cos_theta 0.93     # |cos(theta)| < 0.93
               Vz        10.0     # |Vz| < 10 cm
               Vr        10.0     # Vr < 1 cm (10 mm)
               nChrp     "==2"    # two positive tracks
               nChrn     "==2"    # two negative tracks
               nNet      "==0"    # net charge zero
             }
             .pid(method: :probability) {
               prob_cut 0.001
               identify :proton, against: [:kaon, :pion]   # p+ and anti-p-
               identify :pion,   against: [:kaon, :proton] # pi+ and pi-
               nprp "==1"         # exactly one proton
               nprm "==1"         # exactly one anti-proton
               npip "==1"         # exactly one pi+
               npim "==1"         # exactly one pi-
             }
             .kinematic_fit([:prp, :prm, :pip, :pim]) {
               nominal
               constrain_four_momentum
               chi2_cut 200
             }

alg_pipippbar.with_decay_card(decay_card_pipippbar).apply(sel_pipippbar)
alg_pipippbar.execute_on(data_points + incMC_points + exMC_pipippbar)