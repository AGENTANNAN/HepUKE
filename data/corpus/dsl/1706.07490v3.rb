# =============================================================================
# 1706.07490v3
# Measurement of the cross sections of e+e- -> phi phi omega and
# e+e- -> phi phi phi at c.m. energies from 4.008 to 4.600 GeV (BESIII)
#
# BOSS-side spec: dataset preparation + event selection up to the final
# reconstruction step. The Born cross-section extraction, the ratio
# r_cs = sigma(phi phi omega)/sigma(phi phi phi) and the unbinned fits to the
# recoil-mass spectrum are ROOT-level and are not part of this spec.
#
# Analysis strategy: partial reconstruction. Two phi mesons are reconstructed in
# the prominent phi -> K+ K- mode and the remaining omega (or the third phi) is
# identified through the mass recoiling against the reconstructed phi phi
# system, RM(phi phi). Only part of the final state is detected, so the
# selection uses partial_miss, which replaces the kinematic fit entirely.
#
# Generator: the paper models ISR with KKMC and generates both signals flat in
# phase space (PHSP), so the standard KKMC + psi(4260) top-mother convention is
# used (not ConExc, which additionally has no phi phi phi / phi phi omega mode).
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets — six energy points between sqrt(s) = 4.008 and 4.600 GeV
### Sample names are taken from the BESIII dataset table so that the tabulated
### luminosities match Table 1 of the paper.
### ---------------------------------------------------------------------------
data_4008 = DatasetManager.real_data.find("703_4009")   #  482.0 pb^-1
data_4226 = DatasetManager.real_data.find("703_4230")   # 1091.7 pb^-1
data_4258 = DatasetManager.real_data.find("703_4260")   #  825.7 pb^-1
data_4358 = DatasetManager.real_data.find("703_4360")   #  539.8 pb^-1
data_4416 = DatasetManager.real_data.find("703_4420")   # 1073.6 pb^-1
data_4600 = DatasetManager.real_data.find("703_4600")   #  566.9 pb^-1

data_points = [data_4008, data_4226, data_4258, data_4358, data_4416, data_4600]

# Corresponding inclusive MC samples, for the points available in the table.
incMC_points = [DatasetManager.inclusive_mc.find("703_4009"),
                DatasetManager.inclusive_mc.find("703_4230"),
                DatasetManager.inclusive_mc.find("703_4260"),
                DatasetManager.inclusive_mc.find("703_4360"),
                DatasetManager.inclusive_mc.find("703_4420"),
                DatasetManager.inclusive_mc.find("703_4600")]

### ---------------------------------------------------------------------------
### Decay cards — signal
### ---------------------------------------------------------------------------
# Signal I: e+e- -> phi phi omega.
# The two phi mesons that are reconstructed are forced to K+ K-; the omega is
# never reconstructed (it is identified only through the recoil mass), so it is
# left undecayed — the partial-reconstruction idiom for the untagged particle.
# rec_id_list: 0 psi(4260), 1-2 phi, 3 omega, 4-7 the four kaons.
decay_card_phiphiomega = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi phi omega PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# Signal II: e+e- -> phi phi phi.
# Two phi mesons are reconstructed in K+ K-; the third phi is the untagged
# recoil and is left undecayed. The Alias keeps the untagged phi distinct from
# the two reconstructed ones so that it is not added to the reconstruction list.
# rec_id_list: 0 psi(4260), 1-2 phi, 3 phi_recoil, 4-7 the four kaons.
decay_card_phiphiphi = <<~DECAYCARD
  Alias phi_recoil phi

  Decay psi(4260)
  1.000 phi phi phi_recoil PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Decay cards — peaking backgrounds (estimated with dedicated MC samples)
### ---------------------------------------------------------------------------
# Dominant peaking backgrounds of e+e- -> phi phi phi: when the directly
# produced K+ K- (K+ K- K+ K-) pair is mis-reconstructed as a phi (phi phi),
# e+e- -> K+ K- phi phi and e+e- -> K+ K- K+ K- phi fake the signal. The
# contamination is ~1.0% and 0.1%, respectively.
decay_card_bkg_KKphiphi = <<~DECAYCARD
  Decay psi(4260)
  1.000 K+ K- phi phi PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

decay_card_bkg_KKKKphi = <<~DECAYCARD
  Decay psi(4260)
  1.000 K+ K- K+ K- phi PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# Dominant peaking backgrounds of e+e- -> phi phi omega: the analogous
# e+e- -> K+ K- phi omega and e+e- -> K+ K- K+ K- omega processes.
decay_card_bkg_KKphiomega = <<~DECAYCARD
  Decay psi(4260)
  1.000 K+ K- phi omega PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

decay_card_bkg_KKKKomega = <<~DECAYCARD
  Decay psi(4260)
  1.000 K+ K- K+ K- omega PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC — 0.5 million events per signal process per energy point
### ---------------------------------------------------------------------------
exMCs_phiphiomega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ee_to_phiphiomega_exclusive_mc"   # auto-suffixed per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_phiphiomega
  config.cross_section = :default
end

exMCs_phiphiphi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ee_to_phiphiphi_exclusive_mc"     # auto-suffixed per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_phiphiphi
  config.cross_section = :default
end

# Peaking-background samples (one per energy point each)
exMCs_bkg_KKphiphi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_KKphiphi_exclusive_mc"        # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_KKphiphi
  config.cross_section = :default
end

exMCs_bkg_KKKKphi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_KKKKphi_exclusive_mc"         # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_KKKKphi
  config.cross_section = :default
end

exMCs_bkg_KKphiomega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_KKphiomega_exclusive_mc"      # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_KKphiomega
  config.cross_section = :default
end

exMCs_bkg_KKKKomega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_KKKKomega_exclusive_mc"       # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_KKKKomega
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — e+e- -> phi phi omega
### ---------------------------------------------------------------------------
alg_name_omega = "PhiPhiOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
# Multi-energy cross-section measurement: ECMS is injected per job at run time,
# so it is deliberately NOT set via set_constant here.

selection_omega = Selection.new
  # Four charged tracks: the two phi mesons, each decaying to K+ K-.
  .select_track do
    cos_theta 0.93   # |cos(theta)| < 0.93, polar angle in the MDC
    Vz        10.0   # PCA within +-10 cm along the beam direction
    Vr        1.0    # PCA within 1 cm in the plane perpendicular to the beam
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  # Kaon identification from dE/dx and TOF: L(K) > L(pi).
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==2"
    nkm "==2"
  end
  # Partial reconstruction: both phi mesons are built from their K+ K- daughters,
  # and the omega is undetected, inferred from the recoil against the phi phi
  # system. recID 3 = omega is the missed particle; recIDs 1, 2 and 4-7 (the two
  # phi mesons and their four kaon daughters) are reconstructed.
  # The two best_combination_by_mass calls select the phi phi pairing that
  # minimises sqrt((M_phi1 - M_phi)^2 + (M_phi2 - M_phi)^2), i.e. the minimum
  # Delta M criterion, when more than two phi candidates share no tracks.
  # No recoil-mass window is applied: RM(phi phi) is the fit observable.
  .partial_miss([3]) do
    best_combination_by_mass :phi_psi_4260_1, 1.019461   # PDG phi mass, GeV/c^2
    best_combination_by_mass :phi_psi_4260_2, 1.019461
  end

alg_omega
  .note(:phi_mass_window, "each phi candidate is formed from a pair of oppositely charged
    identified kaons and required to satisfy 1.01 < M(K+K-) < 1.03 GeV/c^2; at least two phi
    candidates with no shared tracks are required per event. When more than two phi candidates
    exist, only the phi phi combination with the minimum
    Delta M = sqrt((M_phi1 - M_phi)^2 + (M_phi2 - M_phi)^2) is kept (M_phi from the PDG) and
    the two phi are randomly labelled phi_1 / phi_2. The 1.01-1.03 GeV/c^2 window itself has no
    dedicated BOSS equivalent outside a fit block and is applied at the ROOT level together
    with the 2D fit of M_phi1 vs M_phi2.")
  .note(:recoil_mass_definition, "the mass recoiling against the reconstructed phi phi system is
    RM(phi phi) = sqrt((E_cm - E_phi phi)^2 - p_phi phi^2), where E_cm is obtained from the
    di-muon process e+e- -> gamma_ISR/FSR mu+ mu- with a precision of 0.02%; the omega signal
    appears as a peak at the omega mass. This variable is the observable of the unbinned
    maximum-likelihood fit and belongs to the ROOT analysis.")
  .note(:recoil_mass_shift_correction, "the measured omega and phi masses in RM(phi phi) are
    shifted left by ~4.5 MeV with respect to the PDG values. The shift is attributed to ISR,
    kaon energy loss, FSR and the E_cm uncertainty and is corrected as an overall shift of
    E_cm, Delta E_cm, calibrated with the control process e+e- -> phi K+ K- (one phi plus one
    charged kaon partially reconstructed); each event is corrected by
    Delta RM(phi phi) = (E_cm - E_phi phi)/RM(phi phi) * Delta E_cm. After the correction the
    fitted omega and phi masses agree with the PDG.")
  .note(:sideband_background, "background without a phi phi pair is estimated from the 2D
    sidebands 0.99 < M(K+K-) < 1.00 GeV/c^2 and 1.04 < M(K+K-) < 1.06 GeV/c^2, as the weighted
    sum of the horizontal and vertical sideband events with the diagonal sidebands subtracted
    to remove double counting; the 2D probability density functions are built from the product
    of two 1D functions (MC-derived phi peak convoluted with a Gaussian, second-order
    polynomial for the non-phi component). Signal leakage into the sidebands is 3%-5%. The
    peaking background from e+e- -> K+ K- phi omega and e+e- -> K+ K+ K- K- omega gives a 1.0%
    systematic uncertainty.")
  .note(:fit_model, "an unbinned maximum-likelihood fit is performed on the corrected
    RM(phi phi) distribution: the signal is modelled by the MC-derived shape and the background
    by a third-order Chebyshev polynomial with parameters fixed from a fit to all samples
    combined. The signal-shape systematic (MC shape convoluted with a free Gaussian) is
    negligible; the background-shape systematic is estimated with a fourth-order polynomial.
    The statistical significance is obtained from the change in likelihood with and without
    the omega component.")
  .note(:efficiency_curve, "the reconstruction efficiency of the two phi mesons depends on their
    production angles, so the PHSP signal MC is re-weighted according to the 2D distribution of
    cos(theta_1) vs cos(theta_2) measured in data; the difference with and without re-weighting
    is taken as the simulation-model systematic uncertainty (2.6%-6.6%).")
  .note(:radiative_correction, "the Born cross section is
    sigma^B = N^obs / (L_int * (1 + delta^r) * (1 + delta^v) * epsilon * B^2), where
    (1 + delta^r) is the radiative correction factor computed in QED with the observed
    energy-dependent line shape as input and interpolated linearly, iterated four times to
    convergence, (1 + delta^v) is the vacuum-polarisation factor from a QED calculation
    (accuracy 0.5%), and B = B(phi -> K+ K-).")
  .with_decay_card(decay_card_phiphiomega)
  .apply(selection_omega)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — e+e- -> phi phi phi
### ---------------------------------------------------------------------------
alg_name_phiphiphi = "PhiPhiPhi"
alg_phiphiphi = Algorithm.new(alg_name_phiphiphi)
alg_phiphiphi.set_header(["#{alg_name_phiphiphi}Alg/#{alg_name_phiphiphi}.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.

selection_phiphiphi = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==2"
    nkm "==2"
  end
  # Partial reconstruction: two phi mesons are rebuilt from K+ K- while the third
  # phi (recID 3 = phi_recoil) is left undetected and inferred from the recoil
  # against the phi phi system. recIDs 1, 2 and 4-7 are reconstructed.
  .partial_miss([3]) do
    best_combination_by_mass :phi_psi_4260_1, 1.019461
    best_combination_by_mass :phi_psi_4260_2, 1.019461
  end

alg_phiphiphi
  .note(:phi_mass_window, "each phi candidate is formed from a pair of oppositely charged
    identified kaons with 1.01 < M(K+K-) < 1.03 GeV/c^2, at least two phi candidates sharing no
    tracks being required; when more candidates exist the pair with the minimum
    Delta M = sqrt((M_phi1 - M_phi)^2 + (M_phi2 - M_phi)^2) is kept. The window and the 2D fit
    of M_phi1 vs M_phi2 are ROOT-level steps.")
  .note(:recoil_mass_definition, "the third phi is identified as a peak at the phi mass in
    RM(phi phi) = sqrt((E_cm - E_phi phi)^2 - p_phi phi^2), with E_cm from the di-muon process
    e+e- -> gamma_ISR/FSR mu+ mu-. The unbinned maximum-likelihood fit to RM(phi phi) is a
    ROOT-level step.")
  .note(:recoil_mass_shift_correction, "the ~4.5 MeV left shift of the measured omega/phi masses
    in RM(phi phi) is corrected as an overall E_cm shift Delta E_cm calibrated with the control
    process e+e- -> phi K+ K-.")
  .note(:sideband_background, "backgrounds are estimated from the 2D phi sidebands
    0.99 < M(K+K-) < 1.00 GeV/c^2 and 1.04 < M(K+K-) < 1.06 GeV/c^2 by a weighted sideband sum
    with the diagonal sidebands subtracted, signal leakage being at the 3%-5% level. The
    dominant peaking backgrounds of this channel, e+e- -> K+ K- phi phi and
    e+e- -> K+ K- K+ K- phi, contribute ~1.0% and 0.1%; 1.0% is taken as the systematic
    uncertainty.")
  .note(:fit_model, "the RM(phi phi) distribution is fitted with a third-order Chebyshev
    polynomial background whose parameters are fixed from all samples combined, the signal
    being described by the MC-derived shape; the background-shape systematic is estimated with
    a fourth-order polynomial and the signal-shape one is negligible.")
  .note(:phi_mass_window_systematic, "the systematic uncertainty from the M(K+K-) requirement is
    neglected, because the mean and width of the Gaussian describing the data/MC difference are
    consistent with zero within three times their uncertainties.")
  .note(:efficiency_curve, "the PHSP signal MC is re-weighted according to the 2D distribution of
    cos(theta_1) vs cos(theta_2) of the two reconstructed phi mesons to reproduce the data
    angular distribution; the resulting efficiency difference (3.4%-10.2%) is the
    simulation-model systematic uncertainty.")
  .note(:radiative_correction, "the Born cross section is
    sigma^B = N^obs / (L_int * (1 + delta^r) * (1 + delta^v) * epsilon * B^2), with the
    radiative correction factor (1 + delta^r) from a QED calculation using the observed line
    shape interpolated linearly and iterated four times to convergence, and the
    vacuum-polarisation factor (1 + delta^v) from QED (accuracy 0.5%).")
  .with_decay_card(decay_card_phiphiphi)
  .apply(selection_phiphiphi)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
root_files_omega = alg_omega.execute_on(
  data_points + incMC_points +
  [exMCs_phiphiomega, exMCs_bkg_KKphiomega, exMCs_bkg_KKKKomega].flatten
)

root_files_phiphiphi = alg_phiphiphi.execute_on(
  data_points + incMC_points +
  [exMCs_phiphiphi, exMCs_bkg_KKphiphi, exMCs_bkg_KKKKphi].flatten
)
