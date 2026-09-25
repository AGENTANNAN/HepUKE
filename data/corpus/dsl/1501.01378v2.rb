# =============================================================================
# BESIII measurement of y_CP in D0-D0bar oscillation using quantum correlations
# in e+e- -> D0 D0bar at sqrt(s) = 3.773 GeV  [arXiv:1501.01378v2]
#
# Doubly-tagged D0 D0bar events: one D is reconstructed in a CP eigenstate
# (the ST tag), the partner D is reconstructed in a semileptonic mode
# (K- e+ nu_e or K- mu+ nu_mu), from the tracks and showers the tag did not use.
# This is a D-tag analysis: single tag side + signal side with a missing neutrino
# (ST + missing pattern).  Two algorithms are needed because the signal-side
# charged multiset differs between the electron and the muon channel.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
psi3770_data  = DatasetManager.real_data.find("712_3773")     # 2.92 fb^-1 at 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # ~10x the data luminosity

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Sub-decay blocks shared by the CP-mode signal cards
sub_ks_pipi = <<~DECAYCARD
  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay
DECAYCARD

sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_eta_gg = <<~DECAYCARD
  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_omega_3pi = <<~DECAYCARD
  Decay omega
  1.0000 pi+ pi- pi0  PHSP;
  Enddecay
DECAYCARD

# psi(3770) -> D0 D0bar with a CP eigenstate on one side and a semileptonic
# decay on the other. The CP final states are charge symmetric, so the CP mode
# can sit on either the D0 or the anti-D0; the card is written with the
# semileptonic decay on the D0 (charge conjugate implicitly covered).
def cp_semileptonic_card(cp_daughters, sub_decays, lepton_daughters)
  <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0  PHSP;
    Enddecay

    Decay D0
    1.0000 #{lepton_daughters}  PHOTOS ISGW2;
    Enddecay

    Decay anti-D0
    1.0000 #{cp_daughters}  PHSP;
    Enddecay

    #{sub_decays}
    End
  DECAYCARD
end

# CP eigenstate tag modes used in this analysis (Table 1):
#   CP-even: K+ K-, pi+ pi-, K_S0 pi0 pi0
#   CP-odd : K_S0 pi0, K_S0 omega
cp_modes = {
  "KK"       => ["K+ K-",          ""],
  "PiPi"     => ["pi+ pi-",        ""],
  "KsPi0Pi0" => ["K_S0 pi0 pi0",   sub_ks_pipi + sub_pi0_gg],
  "KsPi0"    => ["K_S0 pi0",       sub_ks_pipi + sub_pi0_gg],
  "KsOmega"  => ["K_S0 omega",     sub_ks_pipi + sub_omega_3pi + sub_pi0_gg],
}

# Semileptonic signal modes (DT): D0 -> K- e+ nu_e and D0 -> K- mu+ nu_mu
semileptonic_modes = {
  "Kenu"  => ["K- e+ nu_e",   :ep],
  "Kmunu" => ["K- mu+ nu_mu", :mup],
}

# One exclusive MC sample per (CP tag mode x semileptonic mode): these are the
# samples used to evaluate the DT detection efficiencies; the corresponding ST
# efficiencies come from D -> CP, Dbar -> anything samples.
exMC_dt = cp_modes.flat_map do |cp_name, (cp_daughters, cp_subs)|
  semileptonic_modes.map do |sl_name, (sl_daughters, _lepton)|
    DatasetManager.create_exclusive_mc do |config|
      config.sample_name     = "dt_#{cp_name}_#{sl_name}"
      config.related_dataset = psi3770_data
      config.events          = 200_000
      config.decay_card      = cp_semileptonic_card(cp_daughters, cp_subs, sl_daughters)
      config.cross_section   = :default
    end
  end
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — TagAnalysis (D-tag: ST CP tag + semileptonic side)
### ---------------------------------------------------------------------------
# All six CP tag modes are declared on a single tag side; both charm charges are
# scanned (charm-less findSTag overload) because the CP final states are
# charge symmetric.
cp_tag_modes = [:D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toKsPi0, :D0toKsOmega]

### --- Algorithm 1: DT with D0 -> K- e+ nu_e on the signal side ---
alg_kenu_name = "D0CPTagKenu"
alg_kenu = TagAnalysis.new(alg_kenu_name)
alg_kenu.set_header(["#{alg_kenu_name}Alg/#{alg_kenu_name}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(cp_semileptonic_card("K+ K-", "", "K- e+ nu_e"))

alg_kenu.tag_side(:D0) do |t|
  t.modes(*cp_tag_modes)     # CP-eigenstate single tag
end

# Signal side: the K- and the e+ that the tag did not use, plus the neutrino.
# Exactly two oppositely-charged tracks are required (the paper rejects any
# event with a different unused-track multiplicity).
alg_kenu.signal_side do |s|
  s.charged(km: 1, ep: 1)    # exactly two signal-side tracks
  s.require_charge 0         # -1 (K-) + 1 (e+) = 0
  s.missing :nu_e            # massless missing neutrino (semileptonic)
end

# 4C fit: tag D + K- + e+ + nu_e = measured CMS four-momentum.  The tag-side
# mass constraint corresponds to the beam-constrained mass M_BC used in the
# paper.  The published yields come from fits to M_BC (ST) and U_miss (DT), not
# from this fit; the constraint is declared because a tag fit requires it.
alg_kenu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

### --- Algorithm 2: DT with D0 -> K- mu+ nu_mu on the signal side ---
alg_kmunu_name = "D0CPTagKmunu"
alg_kmunu = TagAnalysis.new(alg_kmunu_name)
alg_kmunu.set_header(["#{alg_kmunu_name}Alg/#{alg_kmunu_name}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(cp_semileptonic_card("K+ K-", "", "K- mu+ nu_mu"))

alg_kmunu.tag_side(:D0) do |t|
  t.modes(*cp_tag_modes)
end

alg_kmunu.signal_side do |s|
  s.charged(km: 1, mup: 1)   # exactly two signal-side tracks
  s.require_charge 0         # -1 (K-) + 1 (mu+) = 0
  s.missing :nu_mu           # massless missing neutrino (semileptonic)
end

alg_kmunu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

### ---------------------------------------------------------------------------
### Inexpressible BOSS-side procedures (captured for the systematics step)
### ---------------------------------------------------------------------------
[alg_kenu, alg_kmunu].each do |alg|
  alg.note(:tag_side_selection,
           "ST D candidates: charged tracks satisfy |cos(theta)| < 0.93, closest " \
           "approach to the IP within 1 cm transverse to the beam and within 10 cm " \
           "along the beam. Hadron PID combines the MDC dE/dx and the TOF flight time " \
           "into a likelihood L(h); K candidates require L(K) > L(pi) and pi candidates " \
           "L(pi) > L(K). K_S0 daughters are not subjected to this PID and use a looser " \
           "IP requirement (20 cm along the beam, no transverse constraint).")
     .note(:deltaE_requirements,
           "Mode-dependent deltaE windows (Table 2) are applied on the ST D candidates: " \
           "K+K- [-0.020, 0.020]; pi+pi- [-0.030, 0.030]; K_S0 pi0 pi0 [-0.080, 0.045]; " \
           "K_S0 pi0 [-0.070, 0.040]; K_S0 omega [-0.050, 0.030]; K_S0 eta [-0.040, 0.040] GeV. " \
           "The windows differ per mode and are therefore applied on the stored deltaE in " \
           "the ROOT analysis (a single declared tag-side window cannot express them).")
     .note(:candidate_ranking,
           "Only one candidate per mode per event is accepted: when multiple candidates " \
           "are present the one with the smallest |deltaE| is chosen. The correlation " \
           "between deltaE and M_BC is small, so this does not bias the M_BC background shape.")
     .note(:mBC_signal_region,
           "The ST and DT yields are extracted from the signal region " \
           "1.855 < M_BC < 1.875 GeV/c^2. ST yields come from unbinned maximum-likelihood " \
           "fits to M_BC with a MC-derived signal shape convoluted with a bifurcated " \
           "Gaussian and an ARGUS background; the K_S0 omega mode uses a binned " \
           "least-square fit after omega-sideband subtraction. This is the quantity that " \
           "corresponds to the declared tag-side mass constraint in the fit.")
     .note(:ks_reconstruction,
           "K_S0 candidates are built from pairs of oppositely charged tracks with a " \
           "vertex-constrained fit (vertex chi2 < 100), no PID (both tracks assumed " \
           "pions), with 0.487 < M(pi+pi-) < 0.511 GeV/c^2 and a decay vertex separated " \
           "from the IP by more than two standard deviations.")
     .note(:photon_selection,
           "EMC showers separated from any extrapolated charged track by more than 10 " \
           "standard deviations are photon candidates, with a minimum energy of 25 MeV in " \
           "the barrel (|cos(theta)| < 0.80) and 50 MeV in the end cap " \
           "(0.84 < |cos(theta)| < 0.92); showers in the barrel/end-cap gap are excluded " \
           "and the shower time must be within 700 ns of the event start time.")
     .note(:pi0_eta_reconstruction,
           "pi0 and eta candidates are formed from photon pairs, rejecting pairs with both " \
           "photons in the end cap, with 0.115 < M(gamma gamma) < 0.150 GeV/c^2 (pi0) and " \
           "0.505 < M(gamma gamma) < 0.570 GeV/c^2 (eta); the photon pair is then " \
           "kinematically constrained to the nominal meson mass (1C mass constraint).")
     .note(:omega_reconstruction,
           "omega candidates are reconstructed through omega -> pi+ pi- pi0 with the signal " \
           "region 0.7600 < M(pi+ pi- pi0) < 0.8050 GeV/c^2 and sidebands (0.6000, 0.7300) " \
           "and (0.8300, 0.8525) GeV/c^2; the sidebands are scaled with a factor obtained " \
           "from a fit to the M(pi+ pi- pi0) distribution to subtract the peaking " \
           "background from non-omega D -> K_S0 pi+ pi- pi0 decays.")
     .note(:cosmic_bhabha_veto,
           "For the K+K- and pi+pi- tag modes, events containing only the two tag tracks " \
           "must additionally satisfy: at least one EMC shower separated from the tag " \
           "tracks with energy greater than 50 MeV, and the two tag tracks must not both " \
           "be identified as muons or electrons and (if both have valid TOF times) must " \
           "have a TOF time difference less than 5 ns. This veto is event-level and not " \
           "expressible as a per-candidate selection.")
     .note(:signal_side_track_multiplicity,
           "The semileptonic side must contain exactly two oppositely-charged unused " \
           "tracks satisfying the fiducial requirements; events with any other unused " \
           "track multiplicity are rejected.")
     .note(:semileptonic_pid,
           "K e nu selection: the electron candidate must satisfy L'(e) > 0.001 and " \
           "R'(e) = L'(e)/[L'(e)+L'(pi)+L'(K)] > 0.8 (EMC, dE/dx and TOF combined); when " \
           "both tracks pass, the one with the larger R'(e) is taken as the electron and " \
           "the remaining track must satisfy L(K) > L(pi). K mu nu selection: kaon " \
           "candidates require L(K) > L(pi), the track with the larger L(K) is the K and " \
           "the other is the muon; the muon EMC energy deposit must be less than 0.3 GeV " \
           "and R_l'(e) = L'(e)/[L'(e)+L'(mu)+L'(pi)+L'(K)] must be less than 0.8 to " \
           "suppress D -> K e nu. These likelihood-ratio criteria are not expressible as " \
           "DSL PID cuts.")
     .note(:background_veto,
           "The K mu nu channel requires the K mu invariant mass M(K mu) < 1.65 GeV/c^2 to " \
           "reject D -> K pi backgrounds, and the total energy of unmatched EMC showers " \
           "E_extra < 0.2 GeV to suppress D -> K pi pi0 backgrounds.")
     .note(:umiss_extraction,
           "The DT yields are obtained from unbinned maximum-likelihood fits to " \
           "U_miss = E_miss - c|p_miss| (binned least-square fits with omega-sideband " \
           "subtraction for the modes containing an omega), where " \
           "E_miss = E_beam - E_K - E_l and " \
           "p_miss = -[p_K + p_l + p_ST_hat * sqrt(E_beam^2/c^2 - c^2 m_D^2)]; the " \
           "signal peaks at zero. The K e nu background is a first-order polynomial; the " \
           "K mu nu fit has three components (resolution-corrected K pi pi0 from data, " \
           "K e nu modelled by MC with a fixed 3.5% ratio, and a polynomial for the rest).")
     .note(:no_kinematic_fit_in_paper,
           "The published analysis does NOT perform a kinematic fit: the ST and DT yields " \
           "are extracted from fits to M_BC and U_miss. The 4C fit declared here (with the " \
           "tag-side D0 mass constraint corresponding to M_BC) is only the DSL end-point " \
           "and is not part of the published event selection.")
     .note(:cp_purity,
           "The CP purity of the tag modes is verified by searching for same-CP double-tag " \
           "signals in data: the K_S0 pi0 mode is found to be more than 99% pure and the " \
           "dilution of the C-odd initial state is less than 2% at 90% confidence level. " \
           "Using K+K- as a clean CP-even tag and K_S0 pi0 as a clean CP-odd tag, the " \
           "purities of K_S0 pi0 pi0, K_S0 omega and K_S0 eta are estimated to be larger " \
           "than 89.4%, 93.3% and 93.9% respectively.")
     .note(:charge_conjugate_signal,
           "The CP tag final states are charge symmetric, so both charm charges are " \
           "scanned on the tag side and the declared signal multiset (K- e+ / K- mu+) " \
           "represents the D0 -> K- l+ nu_l combination; the charge-conjugate " \
           "K+ e- / K+ mu- combination is handled symmetrically in the extraction of " \
           "y_CP, which is built from the ratio of the CP-even and CP-odd semileptonic " \
           "branching fractions.")
end

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_kenu.apply
alg_kmunu.apply

alg_kenu.execute_on([psi3770_data, psi3770_incMC] + exMC_dt)
alg_kmunu.execute_on([psi3770_data, psi3770_incMC] + exMC_dt)
