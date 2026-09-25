# Search for charged lepton flavor violating decay J/psi -> e tau
# with tau -> pi pi0 nu_tau
# BESIII: ~10e9 J/psi events at 3.097 GeV

# --- Datasets ---
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# --- Decay card for signal MC ---
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 e- tau+ VLL;
    Enddecay
    Decay tau+
    1.000 pi+ pi0 anti-nu_tau TAUHADNU;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_etau_signal"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

# --- Algorithm ---
alg = Algorithm.new("JpsiETau")
alg.set_header(["JpsiETauAlg/JpsiETau.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })

sel = Selection.new
sel.select_track {
      cos_theta 0.8
      Vz 100.0
      Vr 10.0
      nChrp "==1"
      nChrn "==1"
      nNet "==0"
    }
    .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track 10.0
      nGam ">=2"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
      nep ">=1"
      nem ">=1"
      npip ">=1"
      npim ">=1"
    }
    # pi0 reconstruction via Kalman mass-constrained fit
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
    }
    # Kinematic fit with missing neutrino
    .kinematic_fit([:em, :pip, :pi0]) {
      nominal
      constrain_four_momentum
      miss_track_of(:nu_tau)
      chi2_cut 200
    }

alg.with_decay_card(decay_card_signal).apply(sel)
alg.note(:electron_pid, "Electron ID: CL_e/(CL_e+CL_pi) > 0.95, E/p > 0.8.
  Applied at ROOT level via for_each filter on electron candidates.")
  .note(:electron_momentum, "P_e in [1.009, 1.068] GeV/c for monochromatic
  electron from J/psi -> e tau two-body decay; applied at ROOT level.")
  .note(:recoil_mass, "M_e_recoil in [1.742, 1.811] GeV/c^2 to select
  tau mass region; applied at ROOT level.")
  .note(:missing_energy, "E_miss = E_CMS - E_e - E_pi - E_pi0 > 0.43 GeV
  to suppress fully-reconstructed backgrounds; applied at ROOT level.")
  .note(:umiss_signal_region, "Signal region defined in U_miss = E_miss - c|P_miss|.
  U_miss in [-0.081, 0.112] GeV, corresponding to +/-3 sigma of expected resolution.
  Signal yield obtained by counting events in signal region.")
  .note(:continuum_background, "Continuum background estimated from sqrt(s)=3.08 GeV
  (150 pb-1) and sqrt(s)=3.773 GeV (2.93 fb-1) data, normalized assuming 1/s
  cross-section dependence. Radiative Bhabha is dominant continuum background.")
  .note(:jpsi_background, "Dominant J/psi decay backgrounds: pi+pi-pi0, rho pi,
  omega f2(1270), pbar n pi+. Background normalizations from PDG branching fractions
  and exclusive MC efficiencies. lundcharm modeling uncertainty ~16%.")
  .note(:semiblind_analysis, "Semiblind analysis: ~10% of full data sample randomly
  selected for selection optimization and background study. Signal region unblinded
  after fixing all selection criteria.")
  .note(:two_data_samples, "Two data samples: Sample I (2009+2012, 1.31e9 J/psi)
  and Sample II (2018+2019, 8.70e9 J/psi). Different TOF endcap (upgraded in 2015
  with MRPC, 60 ps resolution). Combined upper limit via profile likelihood.")
  .note(:signal_mc_model, "TAUHADNU generator for tau -> pi pi0 nu based on
  conserved vector currents and chiral Lagrangian. Systematic from alternative
  model (TAUVECTORNU + VSS) evaluated as relative efficiency change.")
  .note(:upper_limit, "UL at 90% C.L.: B(J/psi -> e tau) < 7.5e-8.
  Determined via Bayesian method integrating profile likelihood with
  nuisance parameters for efficiency and background.")

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])