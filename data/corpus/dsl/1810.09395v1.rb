### dataset description ###
# Multi-energy scan from 3.90 to 4.60 GeV — representative example shown for psi(4260)
data_4260 = DatasetManager.real_data.find("703_4260")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# Decay card for Mode I: e+e- -> K_S0 K+ pi- pi0
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K_S0 pi- K+ pi0    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: e+e- -> K_S0 K+ pi- eta
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K_S0 pi- K+ eta    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_KsKpipi0"
    config.related_dataset = data_4260
    config.events = 100000
    config.decay_card = decay_card_modeI
    config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_KsKpieta"
    config.related_dataset = data_4260
    config.events = 100000
    config.decay_card = decay_card_modeII
    config.cross_section = :default
end

### Mode I: e+e- -> K_S0 K+ pi- pi0 ###
alg_modeI = Algorithm.new("KsKpiPi0")
alg_modeI.set_header(["KsKpiPi0Alg/KsKpiPi0.h"])
          .set_constant({"ECMS" => [:double, 4.26]})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta  0.93
              Vz         10.0
              Vr         1.0
              nChrp      ">=2"
              nChrn      ">=2"
              nNet       "==0"
          }
          .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              angle_to_track    20.0
              nGam              ">=2"
          }
          .pid(method: :probability) {
              prob_cut 0.001
              identify :kaon, against: [:pion, :proton]
              nkp ">=1"; nkm ">=1"
          }
          .remove([:kp <= :chrgp, :km <= :chrgn])
          .assign({chrgp: :pip, chrgn: :pim})
          .secondary_vertex_fit([:pip, :pim]) {
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
          }
          .kinematic_fit([:K_S0, :kp, :pim, :gamma, :gamma]) {
              nominal
              constrain_four_momentum
              chi2_cut 60
          }

alg_modeI
  .note(:Ks_mass_window, "K_S0 candidate with |M(pi+pi-) - M_K_S0| < 12 MeV/c^2; secondary vertex fit with decay length > 2 sigma; best K_S0 chosen by smallest chi2 of secondary vertex fit; post-4C-fit signal region M(pi+pi-) in (0.488, 0.508) GeV/c^2")
  .note(:pi0_mass_window, "pi0 candidate with M(gamma gamma) in (0.12, 0.15) GeV/c^2 after 4C kinematic fit; best gamma-gamma pair chosen by smallest chi2_4C; sideband regions defined for background subtraction")
  .note(:kaon_pid_highest_probability, "charged tracks identified as kaon/pion by selecting hypothesis with highest probability from combined TOF and dE/dx; one K+ one K-; remaining positive tracks assigned as pi+, negative as pi-")
  .note(:multi_energy_scan, "cross sections measured at multiple energy points from 3.90 to 4.60 GeV using 5.2 fb^-1 data; ISR correction factor determined iteratively with KKMC; vacuum polarization correction applied; Born cross section computed; upper limits set for Y(4260) and Z_c(3900) decays")
  .note(:sideband_method, "signal yield N_sig = N_A - sum(N_B)/2 + sum(N_C)/4 where A=signal region, B and C = sideband regions; K_S0 sideband (0.463-0.483) U (0.513-0.533) GeV/c^2; pi0 sideband (0.08-0.11) U (0.16-0.19) GeV/c^2")
  .note(:mc_mixing, "data-driven mixing MC sample with intermediate resonances (rho(770), K*(892)) weighted according to momentum distributions in data; used for detection efficiency estimation")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

### Mode II: e+e- -> K_S0 K+ pi- eta ###
alg_modeII = Algorithm.new("KsKpiEta")
alg_modeII.set_header(["KsKpiEtaAlg/KsKpiEta.h"])
           .set_constant({"ECMS" => [:double, 4.26]})

sel_modeII = Selection.new
sel_modeII.select_track {
               cos_theta  0.93
               Vz         10.0
               Vr         1.0
               nChrp      ">=2"
               nChrn      ">=2"
               nNet       "==0"
           }
           .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    20.0
               nGam              ">=2"
           }
           .pid(method: :probability) {
               prob_cut 0.001
               identify :kaon, against: [:pion, :proton]
               nkp ">=1"; nkm ">=1"
           }
           .remove([:kp <= :chrgp, :km <= :chrgn])
           .assign({chrgp: :pip, chrgn: :pim})
           .secondary_vertex_fit([:pip, :pim]) {
               build_virtual_particle(:K_S0).by_minimizing_mass_difference
               remove_used_particle_from_candidate_list
           }
           .kinematic_fit([:K_S0, :kp, :pim, :gamma, :gamma]) {
               nominal
               constrain_four_momentum
               chi2_cut 60
           }

alg_modeII
  .note(:Ks_mass_window, "K_S0 candidate with |M(pi+pi-) - M_K_S0| < 12 MeV/c^2; secondary vertex fit with decay length > 2 sigma; best K_S0 chosen by smallest chi2 of secondary vertex fit; post-4C-fit signal region M(pi+pi-) in (0.488, 0.508) GeV/c^2")
  .note(:eta_mass_window, "eta candidate with M(gamma gamma) in (0.52, 0.58) GeV/c^2 after 4C kinematic fit; best gamma-gamma pair chosen by smallest chi2_4C; sideband regions (0.44-0.50) U (0.60-0.66) GeV/c^2 for background subtraction")
  .note(:kaon_pid_highest_probability, "charged tracks identified as kaon/pion by selecting hypothesis with highest probability from combined TOF and dE/dx; one K+ one K-; remaining positive tracks assigned as pi+, negative as pi-")
  .note(:multi_energy_scan, "cross sections measured at multiple energy points from 3.90 to 4.60 GeV using 5.2 fb^-1 data; ISR correction factor determined iteratively with KKMC; vacuum polarization correction applied; Born cross section computed; upper limits set for Y(4260) and Z_c(3900) decays")
  .note(:sideband_method, "signal yield N_sig = N_A - sum(N_B)/2 + sum(N_C)/4 where A=signal region, B and C = sideband regions; K_S0 sideband (0.463-0.483) U (0.513-0.533) GeV/c^2; eta sideband (0.44-0.50) U (0.60-0.66) GeV/c^2")
  .note(:mc_mixing, "data-driven mixing MC sample with intermediate resonances (rho(770), K*(892)) weighted according to momentum distributions in data; used for detection efficiency estimation")
  .note(:Y4260_Zc3900_search, "upper limits on Y(4260) -> K_S0 K+ pi- pi0/eta and Z_c(3900) -> K_S0 K pi, K_S0 K eta determined by fitting Born cross section line shapes with continuum + BW function")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeI.execute_on([data_4260, incMC_4260, exMC_modeI])
alg_modeII.execute_on([data_4260, incMC_4260, exMC_modeII])