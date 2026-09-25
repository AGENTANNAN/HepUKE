# ============================================================
# Dataset preparation : J/psi -> gamma eta', eta' -> pi+pi-pi0 (Mode I)
#                                      and eta' -> pi0pi0pi0  (Mode II)
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi data @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # matching inclusive MC

# Decay card — Mode I : eta' -> pi+ pi- pi0
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode II : eta' -> pi0 pi0 pi0
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each isospin-violating mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_to_pipipimpi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_to_pi0pi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Mode I : J/psi -> gamma eta', eta' -> pi+ pi- pi0
# ============================================================
alg_name_I = "JpsiGammaEtaPToPiPiPi0"
alg_modeI  = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
            cos_theta 0.93
            Vz        10.0     # |Vz| < 10 cm
            Vr        1.0      # Vr < 1 cm
            nChrp     "==1"    # exactly two oppositely charged tracks
            nChrn     "==1"
            nNet      "==0"    # net charge zero
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025   # 25 MeV in barrel
            energyThreshold_e 0.050   # 50 MeV in endcap
            nGam              ">=3"   # radiative gamma + 2 photons from pi0
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K/p
            npip "==1"
            npim "==1"
          }
         # 1C Kalman fit : build pi0 from a photon pair (chi2 < 25)
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0     ">=1"
          }
         # Nominal 6C fit : 4C energy-momentum + m_pi0 (from Kalman) + m_eta'
         .kinematic_fit([:gamma, :pip, :pim, :pi0]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:etap)
            invariant_mass_of(:gamma, :pi0).out_of(0.733, 0.833)  # veto |M(gamma pi0) - m_omega| < 0.05
            chi2_cut 200   # loose cut; published 25 applied in ROOT
          }
         # Competing hypothesis : 4C fit to gamma gamma gamma pi+ pi- (no pi0)
         .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
            constrain_four_momentum
          }
         # Competing hypothesis : 4C fit to gamma gamma pi+ pi-
         .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
            constrain_four_momentum
          }

alg_modeI.note(:radiative_photon_selection,
               "the radiative photon is taken as the highest-energy photon in the event; the remaining photons are combined into the pi0 candidate")
         .with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Mode II : J/psi -> gamma eta', eta' -> pi0 pi0 pi0
# ============================================================
alg_name_II = "JpsiGammaEtaPToPi0Pi0Pi0"
alg_modeII  = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
             nChrp "==0"   # no charged tracks
             nChrn "==0"
             nNet  "==0"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    10.0
             energyThreshold_b 0.025   # 25 MeV in barrel
             energyThreshold_e 0.050   # 50 MeV in endcap
             nGam              ">=7"   # radiative gamma + 6 photons from 3 pi0
           }
          # 1C Kalman fit : build at least three pi0 from photon pairs (chi2 < 25)
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0     ">=3"
           }
          # 1C Kalman fit : build an eta candidate for the veto / competing fit
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta     ">=1"
           }
          # Nominal 8C fit : 4C energy-momentum + m_eta' + three m_pi0 (from Kalman)
          .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
             invariant_mass_of(:gamma, :pi0).out_of(0.733, 0.833)  # veto |M(gamma pi0) - m_omega| < 0.05
             chi2_cut 200   # loose cut; published 70 applied in ROOT
           }
          # Competing hypothesis : 7C fit to gamma eta pi0 pi0
          .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {
             constrain_four_momentum
           }

alg_modeII.note(:radiative_photon_selection,
                "the radiative photon is taken as the highest-energy photon in the event; the remaining photons are combined into the pi0 candidates")
          .note(:background_veto,
                "gamma-gamma pairs with invariant mass in the eta region [0.52, 0.59] GeV/c2 are vetoed; pi0 candidates with |cos(theta_decay)| > 0.95 are rejected")
          .with_decay_card(decay_card_modeII).apply(sel_modeII)

# ============================================================
# Execute on data, inclusive MC and the two signal MC samples
# ============================================================
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])