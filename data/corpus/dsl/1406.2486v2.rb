# BESIII: e+e- -> p pbar pi0 near the psi(3770) resonance
# Born cross section of e+e- -> p pbar pi0 measured at sqrt(s) = 3.773 GeV, 3.650 GeV and at the
# seven merged points of the psi(3770) line-shape scan (3.736 - 3.813 GeV), and the Born cross
# section of psi(3770) -> p pbar pi0 extracted with continuum-resonance interference.
# Scope: dataset preparation + event selection (BOSS side) only.

### Dataset preparation ###
# Fixed-energy data points
data_3773 = DatasetManager.real_data.find("712_3773")     # 2.9 fb^-1 at the psi(3770) peak
data_3650 = DatasetManager.real_data.find("709_3650")     # 44 pb^-1 at sqrt(s) = 3.650 GeV
# psi(3770) line-shape scan: the 25 small scan sets (3.736 - 3.813 GeV) are merged into 7 sets,
# weighted by luminosity. One representative scan point available as a dataset is used here.
data_scan = DatasetManager.real_data.find("712_3780")     # line-shape scan point, sqrt(s) = 3.780 GeV

incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# Signal decay card: e+e- -> psi(3770) -> p pbar pi0. The intermediate nucleon resonances in
# psi(3770) -> p pbar pi0 are unknown, so the decay is generated with a phase-space model.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive background decay card: psi(3770) -> p pbar pi0 gamma (estimated with 20,000 events).
decay_card_bkg = <<~DECAYCARD
    Decay psi(3770)
    1.000 p+ anti-p- pi0 gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC at each energy point (one sample per point, shared card).
sig_points = [data_3773, data_3650, data_scan]

exMC_signal = DatasetManager.create_exclusive_mc_for(sig_points) do |config|
  config.sample_name   = "pppi0_exclusive_mc"        # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

exMC_bkg = DatasetManager.create_exclusive_mc_for(sig_points) do |config|
  config.sample_name   = "pppi0gamma_exclusive_mc"
  config.events        = 20_000
  config.decay_card    = decay_card_bkg
  config.cross_section = :default
end

### Event selection (BOSS) — e+e- -> p pbar pi0, pi0 -> gamma gamma ###
alg_name = "PpbarPi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
# Multi-energy measurement: ECMS is NOT set here; the CM energy is injected per dataset at run time.

sel = Selection.new
sel.select_track do                  # two charged tracks with net charge zero
      cos_theta 0.93                 # |cos(theta)| < 0.93
      Vz        10.0                 # |Vz| < 10 cm (point of closest approach along the beam)
      Vr        1.0                  # Vr < 1 cm (plane perpendicular to the beam axis)
      nChrp     "==1"                # one positive track  (proton)
      nChrn     "==1"                # one negative track (antiproton)
      nNet      "==0"                # net charge zero
    end
    .select_photon do                # photon candidates for pi0 -> gamma gamma
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025        # E > 25 MeV in the barrel EMC (|cos(theta)| < 0.8)
      energyThreshold_e 0.050        # E > 50 MeV in the endcap EMC (0.86 < |cos(theta)| < 0.92)
      angle_to_track    10.0         # photon must be more than 10 deg from a proton track
      nGam              ">=2"        # at least two photon candidates
    end
    .pid(method: :probability) do    # dE/dx + TOF combined confidence levels
      prob_cut 0.001
      identify :proton, against: [:pion, :kaon]   # CL_p > CL_pi and CL_p > CL_K
      nprp ">=1"                     # at least one proton candidate
      nprm ">=1"                     # at least one antiproton candidate
    end
    # Transverse momentum of the (anti)proton must exceed 300 MeV/c, because the data/MC
    # detection-efficiency difference is large at small transverse momenta.
    .for_each(:prp) do
      define(:pt2) { px * px + py * py }
      where { pt2 < 0.09 }           # pT^2 < (0.3 GeV/c)^2
      remove
    end
    .for_each(:prm) do
      define(:pt2) { px * px + py * py }
      where { pt2 < 0.09 }
      remove
    end
    .select_isolated_photon do       # shower suppression against the (anti)proton tracks
      angle_to_prp_track 10.0        # angle between the photon and the proton > 10 deg
      angle_to_prm_track 30.0        # angle between the photon and the antiproton > 30 deg
      nGam               ">=2"
    end
    # 5C kinematic fit: 4-momentum conservation (4C) plus the pi0 mass constraint on the two
    # photons (1C). When more than two photons are present, all p pbar gamma gamma combinations
    # are iterated and the combination with the smallest chi^2 is selected.
    .kinematic_fit([:prp, :prm, :gamma, :gamma]) do
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 50                    # chi^2 of the 5C fit < 50
    end

alg.note(:pi0_mass_window, "The pi0 candidates are required to lie within a +-3 sigma region of the nominal pi0 mass to further reduce background; this is applied to the photon-pair invariant mass and is characterised offline on the NTuple.")
   .note(:dalitz_efficiency, "Since the intermediate products of psi(3770) -> p pbar pi0 (nucleon resonances) are unknown, the detection efficiency is determined as a function of the Dalitz-plot position using simulated events. The measured cos(theta) distribution of the pi0, which deviates from phase space, is fitted from data and used as input to the signal MC generation.")
   .note(:isr_and_radiative_correction, "Initial state radiation is not included in the efficiency determination: the ISR photons (up to 9% of the beam energy) are accounted for afterwards through the radiative correction factor 1+delta, which also includes vertex corrections and the e+e- self-energy and vacuum-polarisation terms. The observed cross section is converted with sigma_0 = sigma_obs/(1+delta), and 1+delta is iterated together with the cross-section lineshape fit until it is stable.")
   .note(:background_subtraction, "The background from radiative return to the lower-lying J^{PC} = 1^{--} psi(3686) and J/psi resonances (estimated with inclusive MC to be below 0.5%) is subtracted from the number of selected events. Background from D Dbar decays at 3.773 GeV is about 0.015% and is neglected. The data at 3.650 GeV contain a psi(3686) tail contribution of 0.136 +- 0.012 nb, which is taken into account in the cross-section calculation. Other psi(3770) decay channels were estimated to be below 0.4% of the selected events at each energy point and are treated as a systematic uncertainty.")
   .note(:continuum_interference_fit, "The Born cross sections at the nine energy points are fitted with sigma(s) = |sqrt(sigma_con) + sqrt(sigma_psi) * m*Gamma/(s - m^2 + i*m*Gamma) * exp(i*phi)|^2, where sigma_con = C/s^lambda describes the continuum amplitude and phi the relative phase between the resonant and continuum amplitudes. The mass and width of the psi(3770) are fixed to their world-average values. This fit is performed offline on the measured cross sections.")
   .note(:scan_data_merging, "The 25 small data sets of the psi(3770) line-shape scan (3.736 - 3.813 GeV) with varying luminosities are merged into 7 data sets; the CM energy of each merged set is the luminosity-weighted average of its constituent sets, and the error introduced by the merging is treated as a systematic uncertainty.")

alg.with_decay_card(decay_card_signal).apply(sel)

### Execution ###
alg.execute_on([data_3773, data_3650, data_scan,
                incMC_3773, incMC_3650,
                exMC_signal, exMC_bkg].flatten)
