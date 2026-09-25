# Measurement of the D -> K- pi+ strong phase difference in psi(3770) -> D0 D0bar
#   [arXiv:1404.4691]
#
# CP-tagging analysis at sqrt(s) = 3.773 GeV with 2.92 fb^-1. The event is tagged
# by a D0 reconstructed in a CP eigenstate (single tag, ST); the D -> K- pi+
# flavour mode is then reconstructed from the tracks the tag did not use
# (double tag, DT). This is a D-tag analysis: two tag sides of the same species.

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # 2.92 fb^-1 at 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # ~10x the data luminosity

# ---------------------------------------------------------------- decay cards
# Sub-decay blocks shared between the CP-mode signal cards
sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_ks_pipi = <<~DECAYCARD
  Decay K_S0
  1.0000 pi+ pi-  PHSP;
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

sub_rho_pipi = <<~DECAYCARD
  Decay rho0
  1.0000 pi+ pi-  PHSP;
  Enddecay
DECAYCARD

# Build a psi(3770) -> D0 D0bar card with a CP-eigenstate mode on one side and
# the flavour mode D0 -> K- pi+ on the other side.
def dt_signal_card(cp_daughters, sub_decays)
  <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0  PHSP;
    Enddecay

    Decay D0
    1.0000 #{cp_daughters}  PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi-  PHSP;
    Enddecay

    #{sub_decays}
    End
  DECAYCARD
end

# CP eigenstate modes used for tagging (Table 1 of the paper):
#   CP-even (S+): K+ K-, pi+ pi-, K_S0 pi0 pi0, pi0 pi0, rho0 pi0
#   CP-odd  (S-): K_S0 pi0, K_S0 eta, K_S0 omega
dt_cp_modes = {
  "D0toKK"       => ["K+ K-",        ""],
  "D0toPiPi"     => ["pi+ pi-",      ""],
  "D0toKsPi0Pi0" => ["K_S0 pi0 pi0", sub_ks_pipi + sub_pi0_gg],
  "D0toPi0Pi0"   => ["pi0 pi0",      sub_pi0_gg],
  "D0toRhoPi0"   => ["rho0 pi0",     sub_rho_pipi + sub_pi0_gg],
  "D0toKsPi0"    => ["K_S0 pi0",     sub_ks_pipi + sub_pi0_gg],
  "D0toKsEta"    => ["K_S0 eta",     sub_ks_pipi + sub_eta_gg],
  "D0toKsOmega"  => ["K_S0 omega",   sub_ks_pipi + sub_omega_3pi + sub_pi0_gg],
}

# DT signal MC, one sample per CP mode: D -> S+/- on one side, Dbar -> K pi on the
# other. These are the samples used to evaluate the DT detection efficiencies;
# the corresponding ST efficiencies come from D -> S+/- , Dbar -> anything samples.
exMC_dt = dt_cp_modes.map do |mode_name, (daughters, subs)|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "dt_#{mode_name}_KPi"
    config.related_dataset = psi3770_data
    config.events          = 200_000
    config.decay_card      = dt_signal_card(daughters, subs)
    config.cross_section   = :default
  end
end

# Flavour decay card: both sides in the K- pi+ / K+ pi- mode, used to study the
# D0 -> K- pi+ (right-sign) and D0bar -> K- pi+ (wrong-sign) reconstruction.
decay_card_kpi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0  PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+  PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

exMC_kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "dt_KPi_KPi"
  config.related_dataset = psi3770_data
  config.events          = 200_000
  config.decay_card      = decay_card_kpi
  config.cross_section   = :default
end

### Event selection (BOSS) — TagAnalysis (D-tag) ###
alg_name = "D0CPTagKPiStrongPhase"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(dt_signal_card("K_S0 pi0", sub_ks_pipi + sub_pi0_gg))

# Tag side 1 — the CP-eigenstate D. Both charm charges are scanned (charm-less
# findSTag overload), since either D0 or D0bar can carry the CP tag.
alg.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPi0Pi0, :D0toRhoPi0,
          :D0toKsPi0, :D0toKsEta, :D0toKsOmega
end

# Tag side 2 — the flavour D -> K- pi+. Two tag sides of the same species give
# the double-tag (DT) pattern: findDTag pairs the modes and both charges are
# scanned so that K- pi+ and K+ pi- are both covered.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.rank_by :inv
end

# No signal side: in the DT both D mesons are fully accounted for by the tags.

# Kinematic fit over the derived participants (tag1 + tag2). The published
# analysis extracts the yields from M_BC / deltaE rather than from a kinematic
# fit; the constraint is declared here because a tag fit requires it.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)  # tag-side mass constraint
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.note(:deltaE_requirements,
         "Mode-dependent deltaE windows (Table 2) are applied in the analysis: " \
         "K+K- [-0.025, 0.025]; pi+pi- [-0.030, 0.030]; K_S0 pi0 pi0 [-0.080, 0.045]; " \
         "pi0 pi0 [-0.080, 0.040]; rho0 pi0 [-0.070, 0.040]; K_S0 pi0 [-0.070, 0.040]; " \
         "K_S0 eta [-0.040, 0.040]; K_S0 omega [-0.050, 0.030]; K- pi+ [-0.030, 0.030] GeV. " \
         "The boundaries are set at about +/-3 sigma, except for the modes containing a " \
         "pi0, which use (-4 sigma, +3.5 sigma) because of their asymmetric distributions. " \
         "The windows are mode-dependent and therefore not expressible as a single " \
         "tag-side window; they are applied on the stored deltaE in the ROOT analysis.")
   .note(:candidate_ranking,
         "In each event only the combination of D candidates with the least |deltaE| is " \
         "kept per mode for the ST selection, and likewise for the D -> K- pi+ candidate " \
         "in the DT selection.")
   .note(:pi0_eta_reconstruction,
         "pi0 and eta candidates are built from photon pairs requiring at least one photon " \
         "reconstructed in the barrel, with mass windows " \
         "0.115 < m(gamma gamma) < 0.150 GeV/c^2 (pi0) and " \
         "0.505 < m(gamma gamma) < 0.570 GeV/c^2 (eta). The invariant mass of each pair is " \
         "then constrained to the nominal pi0 / eta mass and the four-momentum of the " \
         "candidate is updated according to the fit result (1C mass-constrained fit).")
   .note(:ks_reconstruction,
         "K_S0 candidates are reconstructed from pi+ pi- with a vertex-constrained fit on " \
         "all pairs of oppositely charged tracks and no particle identification. These " \
         "tracks use a looser IP requirement (closest approach within 20 cm along the beam " \
         "direction, no requirement in the transverse plane). The vertex fit chi2 must be " \
         "less than 100. A second fit constrains the K_S0 momentum to point back to the IP " \
         "and the resulting flight length must satisfy L/sigma_L > 2. Finally the pi+ pi- " \
         "invariant mass must lie in (0.487, 0.511) GeV/c^2 (three times the experimental " \
         "mass resolution).")
   .note(:cosmic_bhabha_veto,
         "In the K+K- and pi+pi- modes the cosmic-ray and Bhabha backgrounds are removed " \
         "by requiring (a) the two CP-tag tracks to have a TOF time difference less than " \
         "5 ns and to be inconsistent with a muon pair or an electron-positron pair, and " \
         "(b) at least one EMC shower (other than those from the CP-tag tracks) with " \
         "energy greater than 50 MeV, or at least one additional charged track in the MDC.")
   .note(:resonance_mass_windows,
         "The rho0 pi0 and K_S0 omega modes require 0.60 < m(pi+ pi-) < 0.95 GeV/c^2 and " \
         "0.72 < m(pi+ pi- pi0) < 0.84 GeV/c^2 to identify the rho and omega candidates. " \
         "In the K_S0 pi0 mode the D0 -> rho pi background is negligible once the " \
         "L/sigma_L > 2 decay-length requirement is applied.")
   .note(:st_yield_fit,
         "The ST yields n_S+/- are extracted from maximum-likelihood fits to the M_BC " \
         "distributions, with the signal modelled by the reconstructed MC signal shape " \
         "convoluted with a smearing Gaussian and the background by an ARGUS function.")
   .note(:dt_yield_fit,
         "The DT yields n_(K pi, S+/-) are extracted from two-dimensional maximum-likelihood " \
         "fits to M_BC(S+/-) vs M_BC(K pi). The signal shapes are taken from MC simulation; " \
         "the background shapes contain the continuum background and a mis-partitioning " \
         "background in which final-state particles are interchanged between the two D " \
         "candidates.")
   .note(:cp_purity,
         "The CP purity f_S of each ST mode is determined: 98.5% for K_S0 pi0 and almost " \
         "100% for K_S0 eta, from the K_S0 mass sidebands [0.470, 0.477] and " \
         "[0.521, 0.528] GeV/c^2. For K_S0 omega, K_S0 pi0 pi0 and rho0 pi0 the purity is " \
         "measured directly from data using same-CP double-tag combinations (S', S) with " \
         "the clean CP tags K_S0 pi0 (S'-) and K+K- (S'+).")
   .note(:efficiency_ratio_systematics,
         "The efficiencies of the CP-tagged and K pi-tagged sides largely cancel in the " \
         "asymmetry. The residual data-MC difference of the double-tag to single-tag yield " \
         "ratio, Delta_S+/- = Delta(eps_S+/- / eps_(K pi, S+/-)), is studied with control " \
         "samples and found to be at the 1% level per CP-tag mode, giving a systematic " \
         "uncertainty of 0.2 x 10^-2 on A_CP(K pi).")
   .note(:external_inputs,
         "A_CP(K pi) is converted into cos(delta_K pi) using the external inputs " \
         "r^2 = (3.50 +/- 0.04) x 10^-3 and y = (6.7 +/- 0.9) x 10^-3 from HFAG and " \
         "R_WS = (3.80 +/- 0.05) x 10^-3 from the PDG.")
   .note(:st_inclusive_mc,
         "In addition to the DT signal samples, the single-tag efficiencies are evaluated " \
         "with MC samples of D -> S+/- with the opposite D decaying inclusively " \
         "(Dbar -> X); such an inclusive decay of the partner D is not expressible in an " \
         "EvtGen decay card and is produced outside this spec.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_kpi] + exMC_dt)
