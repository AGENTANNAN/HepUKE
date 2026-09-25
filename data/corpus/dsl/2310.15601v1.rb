# Paper: 2310.15601v1
# Title: Study of doubly Cabibbo-suppressed Ds+ → K+ K+ π- and Ds+ → K+ K+ π- π0
# Analysis type: TAG-BASED (Ds DT)
# Datasets: BOSS 703 + 705 scan (8 energy points, 4.128-4.226 GeV)

# === DATASET PREPARATION ===

# Build the dataset list for 8 CM energy points
# 4.128, 4.157 GeV → BOSS 705
# 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV → BOSS 703
data_4128 = DatasetManager.real_data.find("705_4130")
data_4157 = DatasetManager.real_data.find("705_4160")
data_4178 = DatasetManager.real_data.find("703_4180")
data_4189 = DatasetManager.real_data.find("703_4190")
data_4199 = DatasetManager.real_data.find("703_4200")
data_4209 = DatasetManager.real_data.find("703_4210")
data_4219 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

all_data = [data_4128, data_4157, data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]

# Inclusive MC
inc_mc_4130 = DatasetManager.load_inclusive_mc.find("705_4130")
inc_mc_4160 = DatasetManager.load_inclusive_mc.find("705_4160")
inc_mc_4180 = DatasetManager.load_inclusive_mc.find("703_4180")
inc_mc_4190 = DatasetManager.load_inclusive_mc.find("703_4190")
inc_mc_4200 = DatasetManager.load_inclusive_mc.find("703_4200")
inc_mc_4210 = DatasetManager.load_inclusive_mc.find("703_4210")
inc_mc_4220 = DatasetManager.load_inclusive_mc.find("703_4220")
inc_mc_4230 = DatasetManager.load_inclusive_mc.find("703_4230")

# Exclusive signal MC: e+e- → Ds*+ Ds- → γ/π0 Ds+ Ds-
# Ds+ → K+ K+ π- (DCS signal mode 1) or Ds+ → K+ K+ π- π0 (DCS signal mode 2)
# Generated at all 8 energy points
exc_mc_sig1 = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.decay_card <<~DECAYCARD
    Decay D_s+
    1.0 K+ K+ pi-  PHSP;
    Enddecay
  DECAYCARD
end

exc_mc_sig2 = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.decay_card <<~DECAYCARD
    Decay D_s+
    1.0 K+ K+ pi- pi0  PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma  PHSP;
    Enddecay
  DECAYCARD
end

# === TAG ANALYSIS ===

alg = TagAnalysis.new("Ds_DCS_DT")
alg.set_header(["Ds_DCS_DTAlg/Ds_DCS_DT.h"])

# ECMS values per energy point (approximate averages)
# Note: The actual ECMS is per-run from the DB via beam_energy :db
# Using representative values for the constant
alg.set_constant({ "ECMS" => [:double, 4.178] })

# --- Tag side 1: Ds- ST reconstruction ---
# Paper Sec IV lists 11 hadronic tag modes. Modes with available DTagAlg channel
# names are included. Modes involving η/η'/ρ sub-decays that lack dedicated
# DTagAlg symbols are noted as unavailable.
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi       # K+ K- π- (CF tag mode)

  # Additional modes with available DTagAlg symbols:
  # Note: DTagAlg channel names map final-state particles, not intermediate
  # resonances. Modes involving η → γγ/π+π-π0, η' sub-decays, and ρ decays
  # are typically covered by :DstoPiEta or similar final-state combinations,
  # but may not exist for all sub-decay selections.
  t.charm -1               # Pin to Ds- (anti-charm)
end

# The paper also uses these tag modes which lack dedicated DTagAlg channel symbols
# in the frozen BOSS 7.0.6/7.1.2 releases:
alg.note(:tag_mode_unavailable,
  "Tag modes without DTagAlg channel symbols: " \
  "K+K-pi-pi0, pi-pi+pi-, KS0K-, KS0K+pi-pi-, " \
  "eta(gamma gamma)pi-, eta(pi+pi-pi0)pi-, " \
  "eta'(pi+pi-eta)pi-, eta'(gamma rho0)pi-, " \
  "eta(gamma gamma)rho-, eta(pi+pi-pi0)rho-")

# --- Tag side 2: Ds+ signal reconstruction (DT) ---
# The signal Ds+ is reconstructed from the remaining tracks/showers
# not used by the tag Ds-. DT method uses the transition γ/π0 from Ds*.
alg.tag_side(:Ds) do |t|
  # Use the same tag modes (DT = reconstruct both sides from same species)
  t.modes :DstoKKPi       # K+ K- π+ (CF, charged conjugate)
  t.charm 1                # Pin to Ds+
  t.rank_by :inv           # DT ranking by invariant mass difference
end

# --- Signal side: tracks not used by either tag ---
# Two signal channels are analyzed:
# (1) Ds+ → K+ K+ π- (requires 3 charged tracks: 2 K+, 1 π-)
# (2) Ds+ → K+ K+ π- π0 (requires 3 charged tracks + π0 → 2γ)
# Both share the same charged content; π0 selection is done at ROOT level
alg.signal_side do |s|
  s.charged(kp: 2, pim: 1)   # 2 K+, 1 π- not used by tags
  s.photons 0..48             # 0 for mode 1, >=2 for mode 2
end

# --- Kinematic fit (6C: 4C + mass constraints on both Ds tags) ---
# The DT fit constrains total 4-momentum and both Ds masses
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200             # loose cut; optimal in ROOT
end

# --- DTag reconstruction defaults ---
# Default: local true, beam_energy :db, trim_mode_lists true
# These are acceptable for this analysis.

# --- Transition photon/π0 selection note ---
# The paper loops over unused γ and π0 candidates under both
# Ds* → γ Ds and Ds* → π0 Ds hypotheses, choosing min |ΔE|.
# This logic is handled by the DTagAlg in BOSS.
alg.note(:transition_selection,
  "Loop over unused gamma and pi0 candidates; " \
  "choose combination with minimum |Delta E| under " \
  "Ds*+ → gamma Ds+ and Ds*+ → pi0 Ds+ hypotheses")

# --- M_BC signal region windows per energy point ---
# Applied in ROOT analysis; stored unconditionally by DTagAlg
alg.note(:mbc_windows,
  "M_BC windows per energy: 4.128[2.010,2.061], 4.157[2.010,2.070], " \
  "4.178[2.010,2.073], 4.189[2.010,2.076], 4.199[2.010,2.079], " \
  "4.209[2.010,2.082], 4.219[2.010,2.085], 4.226[2.010,2.088] GeV")

alg.with_decay_card(<<~DECAYCARD)
  Decay D_s+
  1.0 K+ K+ pi-  PHSP;
  Enddecay
DECAYCARD

alg.apply
alg.execute_on(all_data)