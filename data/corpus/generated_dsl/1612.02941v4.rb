# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive J/psi MC

# ---------------------------------------------------------------------------
# Decay cards (EvtGen format, EvtGen particle names)
# ---------------------------------------------------------------------------

# Signal chain 1: J/psi -> gamma eta_c, eta_c -> phi phi, phi -> K+ K-
decay_card_signal_phiphi = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 phi phi PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Signal chain 2: J/psi -> gamma eta_c, eta_c -> omega phi, omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
decay_card_signal_omegaphi = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 omega phi PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background for the phi-phi mode: J/psi -> gamma phi K+ K-  (phi -> K+ K-)
decay_card_bg_phiphi_phiKK = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma phi K+ K- PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Background for the phi-phi mode: J/psi -> gamma K+ K- K+ K-  (non-resonant four-kaon)
decay_card_bg_phiphi_4K = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma K+ K- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Background for the omega-phi mode: J/psi -> eta' phi, eta' -> gamma omega
decay_card_bg_omegaphi_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 eta' phi PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma omega PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples
# ---------------------------------------------------------------------------

# 500k signal MC: J/psi -> gamma eta_c, eta_c -> phi phi
exMC_sig_phiphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etac_phiphi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal_phiphi
  config.cross_section   = :default
end

# 500k signal MC: J/psi -> gamma eta_c, eta_c -> omega phi
exMC_sig_omegaphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etac_omegaphi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal_omegaphi
  config.cross_section   = :default
end

# 200k background MC for the phi-phi mode: J/psi -> gamma phi K+ K-
exMC_bg_phiphi_phiKK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_phi_KK_bkg"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bg_phiphi_phiKK
  config.cross_section   = :default
end

# 200k background MC for the phi-phi mode: J/psi -> gamma K+ K- K+ K-
exMC_bg_phiphi_4K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_4K_bkg"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bg_phiphi_4K
  config.cross_section   = :default
end

# 200k background MC for the omega-phi mode: J/psi -> eta' phi, eta' -> gamma omega
exMC_bg_omegaphi_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_etap_phi_bkg"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bg_omegaphi_etap
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ===========================================================================
# Mode I: J/psi -> gamma eta_c, eta_c -> phi phi, phi -> K+ K-
# ===========================================================================
alg_name_phiphi = "JpsiGammaEtacPhiPhi"
alg_phiphi = Algorithm.new(alg_name_phiphi)
alg_phiphi.set_header(["#{alg_name_phiphi}Alg/#{alg_name_phiphi}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})   # ECMS = 3.097 GeV

sel_phiphi = Selection.new
  .select_track {                 # Charged track selection
     cos_theta 0.93               # |cos(theta)| < 0.93
     Vz        10.0               # |Vz| < 10 cm
     Vr        1.0                # Vr < 1 cm
     nChrp     "==2"              # exactly two positive tracks
     nChrn     "==2"              # exactly two negative tracks
     nNet      "==0"              # net charge zero
  }
  .select_photon {                # Photon selection
     tdc_emc_start     0          # EMC timing window 0 - 700 ns
     tdc_emc_end       14         # (unit 50 ns -> 14 * 50 ns = 700 ns)
     angle_to_track    10.0       # at least 10 deg from the nearest charged track
     energyThreshold_b 0.025      # > 25 MeV in the EMC barrel
     energyThreshold_e 0.050      # > 50 MeV in the EMC endcap
     nGam              ">=1"      # at least one photon
  }
  .assign({:chrgp => :kp, :chrgn => :km})   # no PID: take the four tracks directly as K+ K- K+ K-
  .kinematic_fit([:gamma, :kp, :kp, :km, :km]) {   # 4C fit to gamma K+ K- K+ K-
     nominal                    # nominal fit: corrected four-momenta are saved
     constrain_four_momentum    # 4C energy-momentum constraint
     chi2_cut 200               # loose chi^2 cut; optimal tight cut applied in ROOT
  }

# The two phi candidates are fixed by minimising the sum of the K+K- mass deviations
# (candidate pairing before/at the fit) - not expressible by a dedicated DSL method.
alg_phiphi.note(:best_combination_selection,
  "the two phi candidates in the phi-phi mode are fixed by minimising the sum of the " \
  "|M(K+K-) - M_phi| deviations over the two K+K- pairings; the nominal 4C kinematic fit " \
  "then keeps the remaining candidates with the smallest chi^2")

alg_phiphi.with_decay_card(decay_card_signal_phiphi).apply(sel_phiphi)

# ===========================================================================
# Mode II: J/psi -> gamma eta_c, eta_c -> omega phi,
#          omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
# ===========================================================================
alg_name_omegaphi = "JpsiGammaEtacOmegaPhi"
alg_omegaphi = Algorithm.new(alg_name_omegaphi)
alg_omegaphi.set_header(["#{alg_name_omegaphi}Alg/#{alg_name_omegaphi}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})  # ECMS = 3.097 GeV

sel_omegaphi = Selection.new
  .select_track {                 # Charged track selection
     cos_theta 0.93               # |cos(theta)| < 0.93
     Vz        10.0               # |Vz| < 10 cm
     Vr        1.0                # Vr < 1 cm
     nChrp     "==2"              # exactly two positive tracks
     nChrn     "==2"              # exactly two negative tracks
     nNet      "==0"              # net charge zero
  }
  .select_photon {                # Photon selection
     tdc_emc_start     0          # EMC timing window 0 - 700 ns
     tdc_emc_end       14         # (unit 50 ns -> 14 * 50 ns = 700 ns)
     angle_to_track    10.0       # at least 10 deg from the nearest charged track
     energyThreshold_b 0.025      # > 25 MeV in the EMC barrel
     energyThreshold_e 0.050      # > 50 MeV in the EMC endcap
     nGam              ">=3"      # at least three photons
  }
  .pid(method: :probability) {    # Probability PID
     prob_cut 0.001               # PID probability > 0.001
     identify :kaon, against: [:pion, :proton]   # one K+ and one K-
     nkp "==1"
     nkm "==1"
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])   # remove the identified kaons from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim}) # remaining positive/negative tracks are pi+ / pi-
  .kalman_kinematic_fit([:gamma, :gamma]) { # 1C mass-constrained fit to reconstruct pi0 from a gamma-gamma pair
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 25                 # chi^2 < 25
     npi0 ">=1"                  # at least one pi0 (best gamma-gamma pair by smallest chi^2)
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {  # 4C fit to gamma K+ K- pi+ pi- pi0
     nominal                     # nominal fit: corrected four-momenta are saved
     constrain_four_momentum     # 4C energy-momentum constraint
     chi2_cut 200                # loose chi^2 cut; optimal tight cut applied in ROOT
     invariant_mass_of(:kp, :km).within(1.0115, 1.0275)      # |M(K+K-) - M_phi| < 0.008 GeV/c^2
     invariant_mass_of(:pip, :pim, :pi0).within(0.7527, 0.8127)  # |M(pi+pi-pi0) - M_omega| < 0.03 GeV/c^2
  }

alg_omegaphi.with_decay_card(decay_card_signal_omegaphi).apply(sel_omegaphi)

### Execution ###
root_files_phiphi = alg_phiphi.execute_on([
  jpsi_data, jpsi_incMC,
  exMC_sig_phiphi, exMC_bg_phiphi_phiKK, exMC_bg_phiphi_4K
])

root_files_omegaphi = alg_omegaphi.execute_on([
  jpsi_data, jpsi_incMC,
  exMC_sig_omegaphi, exMC_bg_omegaphi_etap
])