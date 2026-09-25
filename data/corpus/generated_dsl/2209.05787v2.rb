# ============================================================================
# BOSS / dataset-preparation DSL for the search of the baryon- and
# lepton-number-violating decays D± → (anti-)n e± at √s = 3.773 GeV.
# The event is tagged by a fully reconstructed hadronic D (single tag), and the
# opposite D is searched for in the e± + missing-(anti-)neutron final state.
# Structure = one tag_side + signal_side with a missing declaration  →  "ST + missing".
# ============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data, √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC sample
# (no exclusive-MC sample size was specified, so only data + inclusive MC are used)

# Decay card: ψ(3770) → D+D−, with the Δ|B−L| = 0 and Δ|B−L| = 2 signal decays.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 anti-n0 e+ PHSP;
    Enddecay

    Decay D-
    1.0000 n0 e- PHSP;
    Enddecay

    End
DECAYCARD

### Tag-based event selection (BOSS) ###
alg_name = "DNeTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# ---- single-tag D± : six hadronic tag modes. charm is left unspecified so both
#      the D+ modes and their charge conjugates (the D− modes) are scanned. ----
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  # ST tag-side windows (store-not-cut is lifted only because they were requested)
  t.window :deltaE, min: -0.055, max: 0.040   # mode-dependent ΔE, union window (see note)
  t.window :mBC,    min: 1.863,  max: 1.877   # M_BC window from the ARGUS fit
end

# ---- signal side: exactly one electron, zero photons, (anti-)neutron missing ----
alg.signal_side do |s|
  s.photons 0          # zero signal-side photons
  s.charged(ep: 1)     # exactly one electron, no additional charged tracks
  s.require_charge(1)  # net signal-side charge ±1
  s.missing :n0        # (anti-)neutron treated as a single missing particle
end

# ---- 2C kinematic fit (convergence only) ----
alg.fit do |f|
  f.constrain_four_momentum                                           # four-momentum conservation
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)     # tag mass → m(D+)
  f.invariant_mass_of(:ep, :n0).constrain_to_nominal_mass_of(:Dplus)  # (e + missing) → m(D∓)
  f.chi2_cut 9999                                                     # only require convergence
end

# ---- BOSS-side procedures that cannot be expressed in the DSL ----
alg.note(:delta_e_mode_dependent,
  "The ST ΔE requirement is mode dependent: |ΔE| < 25 MeV for the tag modes without a π0 " \
  "(K−π+π+, K_S0π+, K_S0π+π+π−, K+K−π+) and −55 < ΔE < +40 MeV for the π0 modes " \
  "(K−π+π+π0, K_S0π+π0). Only one ΔE window per tag side can be declared, so the union " \
  "window (−55, +40) MeV is used here and the mode-specific window is re-applied in ROOT.")

alg.note(:signal_charge_conjugate,
  "This spec encodes the D− tag / D+ → (anti-)n e+ configuration (signal electron e+). The " \
  "charge-conjugate configuration (D+ tag, signal e−, D− → n e−) is identical apart from the " \
  "tagged charm (-1 → +1) and charged(em: 1); the single-tag modes and windows are unchanged.")

alg.note(:neutron_mass,
  "The (anti-)neutron is reconstructed as a single missing particle; its mass is treated as " \
  "unknown (not imposed as a fit constraint), so the only mass constraints are tag → m(D) and " \
  "(e + n) → m(D).")

# ---- render the algorithm and run on data + inclusive MC ----
alg.apply
root_files = alg.execute_on([data_3773, incMC_3773])