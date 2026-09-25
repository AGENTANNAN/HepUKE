# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data (BOSS 7.0.8, 3.097 GeV)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC sample

# Decay card for the REFERENCE mode: J/psi -> phi eta, phi -> K+ K-, eta -> gamma gamma
decay_card_kk = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the SIGNAL (lepton-number-violating) mode:
# J/psi -> phi eta, phi -> pi+ pi+ e- e-, eta -> gamma gamma
decay_card_lnv = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta PHSP;
    Enddecay

    Decay phi
    1.0000 pi+ pi+ e- e- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the reference mode: 1,000,000 events
exMC_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_kk_eta_gg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_kk
  config.cross_section   = :default
end

# Exclusive MC for the signal mode: 500,000 events
exMC_lnv = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_pipi_ee_eta_gg"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_lnv
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------
# Reference mode: J/psi -> phi eta, phi -> K+ K-, eta -> gamma gamma
# ------------------------------------------------------------------
alg_ref_name = "PhiEtaKK"
alg_ref = Algorithm.new(alg_ref_name)
alg_ref.set_header(["#{alg_ref_name}Alg/#{alg_ref_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})           # CMS energy 3.097 GeV
       .set_alias({"std::vector<double>" => "Vdouble"})

ref_selection = Selection.new
ref_selection.select_track {                 # charged track selection (common)
                 cos_theta 0.93              # |cos(theta)| < 0.93
                 Vz        10.0              # |Vz| < 10 cm
                 Vr        1.0               # Vr < 1 cm
                 nChrp     "==1"             # exactly two charged tracks ...
                 nChrn     "==1"
                 nNet      "==0"             # ... with net charge zero
               }
               .select_photon {              # photon selection (common)
                 tdc_emc_start     0         # EMC TDC window 0-700 ns
                 tdc_emc_end       14
                 angle_to_track    20.0      # opening angle > 20 deg w.r.t. nearest charged track
                 energyThreshold_b 0.025     # 25 MeV (barrel)
                 energyThreshold_e 0.050     # 50 MeV (endcap)
                 nGam              ">=2"     # at least two photons
               }
               .pid(method: :probability) {  # PID by the probability method
                 prob_cut 0.001              # CL > 0.001
                 identify :kaon, against: [:pion]   # kaon vs pion (both charges)
                 nkp "==1"                   # one K+
                 nkm "==1"                   # one K-
               }
               .kinematic_fit([:kp, :km, :gamma, :gamma]) {   # 4C fit to K+ K- gamma gamma
                 nominal                     # nominal fit
                 constrain_four_momentum     # 4C energy-momentum constraint
                 chi2_cut 200
               }

alg_ref.with_decay_card(decay_card_kk).apply(ref_selection)

# ------------------------------------------------------------------
# Signal mode: J/psi -> phi eta, phi -> pi+ pi+ e- e-, eta -> gamma gamma
# ------------------------------------------------------------------
alg_sig_name = "PhiEtaLNV"
alg_sig = Algorithm.new(alg_sig_name)
alg_sig.set_header(["#{alg_sig_name}Alg/#{alg_sig_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})           # CMS energy 3.097 GeV
       .set_alias({"std::vector<double>" => "Vdouble"})

sig_selection = Selection.new
sig_selection.select_track {                 # charged track selection (common)
                 cos_theta 0.93              # |cos(theta)| < 0.93
                 Vz        10.0              # |Vz| < 10 cm
                 Vr        1.0               # Vr < 1 cm
                 nChrp     "==2"             # exactly four charged tracks ...
                 nChrn     "==2"
                 nNet      "==0"             # ... with net charge zero
               }
               .select_photon {              # photon selection (common)
                 tdc_emc_start     0         # EMC TDC window 0-700 ns
                 tdc_emc_end       14
                 angle_to_track    20.0      # opening angle > 20 deg w.r.t. nearest charged track
                 energyThreshold_b 0.025     # 25 MeV (barrel)
                 energyThreshold_e 0.050     # 50 MeV (endcap)
                 nGam              ">=2"     # at least two photons
               }
               .pid(method: :probability) {
                 prob_cut 0.001              # CL > 0.001
                 # high-momentum tracks (p > 1.0 GeV) treated as leptons;
                 # electron if the EMC energy deposit exceeds 0.6 GeV, otherwise muon
                 identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                treat_as_electron_if_energy_above: 0.6
                 identify :pion, against: [:kaon]   # pion vs kaon
                 npip ">=2"                  # at least two pi+
                 nlp  ">=1"                  # one positive lepton (e+)
                 nlm  ">=1"                  # one negative lepton (e-)
               }
               # nominal 4C fit to pi+ pi+ e- e- gamma gamma
               .kinematic_fit([:pip, :pip, :lm, :lm, :gamma, :gamma]) {
                 nominal
                 constrain_four_momentum
                 # phi window: 0.99 < M(pi+pi+e-e-) < 1.04 GeV
                 invariant_mass_of(:pip, :pip, :lm, :lm).within(0.99, 1.04)
                 # eta window: 0.52 < M(gamma gamma) < 0.57 GeV
                 invariant_mass_of(:gamma, :gamma).within(0.52, 0.57)
                 chi2_cut 30
               }
               # ------------------------------------------------------------------
               # Six competing 4C hypotheses: their chi2 values are stored (no
               # nominal, no chi2 cut) and the "worse than signal" veto is applied
               # later at the ROOT level.
               # ------------------------------------------------------------------
               .assign({:chrgp => :kp, :chrgn => :km})  # reinterpret the tracks as kaons
               .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) { constrain_four_momentum }       # K+K-pi+pi-gammagamma
               .kinematic_fit([:pip, :pip, :pim, :pim, :gamma, :gamma]) { constrain_four_momentum }     # pi+pi+pi-pi-gammagamma
               .kinematic_fit([:kp, :km, :lp, :lm, :gamma, :gamma]) { constrain_four_momentum }         # K+K-e+e-gammagamma
               .kinematic_fit([:kp, :km, :lp, :pim, :gamma, :gamma]) { constrain_four_momentum }        # K+K-e-pi-gammagamma
               .kinematic_fit([:kp, :km, :pip, :lp, :lm, :gamma, :gamma]) { constrain_four_momentum }   # K K pi e+e- gammagamma
               .kinematic_fit([:kp, :km, :pip, :lp, :pim, :gamma, :gamma]) { constrain_four_momentum }  # K K pi e pi gammagamma

# Opening-angle veto against gamma conversions (not expressible with the current DSL:
# the angle_between / cos_theta_between primitives are not implemented yet)
alg_sig.note(:background_veto,
             "reject events with opening angle theta(pi, e) < 8 degrees between the signal pion and " \
             "electron candidate to suppress photon conversions gamma -> e+e-; the angular variable " \
             "is not expressible in the current DSL and is evaluated in the analysis step")

alg_sig.with_decay_card(decay_card_lnv).apply(sig_selection)

### Execution ###
root_files_ref = alg_ref.execute_on([jpsi_data, jpsi_incMC, exMC_kk])
root_files_sig = alg_sig.execute_on([jpsi_data, jpsi_incMC, exMC_lnv])