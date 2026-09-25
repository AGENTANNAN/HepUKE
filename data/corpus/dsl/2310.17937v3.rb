# Paper: 2310.17937v3
# Title: Observation of the Anomalous Shape of X(1840) in J/ψ → γ 3(π+π-)
# Analysis type: ORDINARY
# Dataset: J/ψ (BOSS 708, 3.097 GeV)

# === DATASET PREPARATION ===

data_jpsi = DatasetManager.real_data.find("708_3097")
inc_mc = DatasetManager.load_inclusive_mc.find("708_3097")

# Exclusive MC: J/ψ → γ X(1840) → γ 3(π+π-), phase space
exc_mc = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~DECAYCARD
    Decay J/psi
    1.0 gamma X1840  PHSP;
    Enddecay
    Decay X1840
    1.0 pi+ pi- pi+ pi- pi+ pi-  PHSP;
    Enddecay
  DECAYCARD
  c.related_dataset data_jpsi
end

# === ALGORITHM & SELECTION ===

alg = Algorithm.new("Jpsi_gamma_6pi", version: '00-00-01')
alg.set_header(["Jpsi_gamma_6piAlg/Jpsi_gamma_6pi.h"])
alg.set_constant({ "ECMS" => [:double, 3.097] })

sel = alg.selection

# --- Track selection ---
# Require exactly 6 charged tracks with zero net charge
sel.select_track do |t|
  t.nTot 6
  # All charged tracks assumed to be pions (no specific PID required per paper)
end

# Note: The paper does NOT apply pion PID — all charged tracks are treated as pions.
# The net charge is zero naturally when 6 tracks are balanced as 3π+ and 3π-.

# --- Photon selection ---
sel.select_photon do |ph|
  ph.nMin 1           # at least 1 photon
end

# --- 4C kinematic fit under J/ψ → γ 3(π+π-) hypothesis ---
# The fit constrains total 4-momentum to ECMS.
# The paper requires χ²_4C < 30 and the best combination (min χ²_4C) retained.
# Post-fit: χ²(J/ψ→γ 6π) < χ²(J/ψ→γγ 6π) to suppress γγ 6π background — applied in ROOT.
sel.kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pip, :pim]) do |fit|
  fit.constrain_four_momentum
  fit.chi2_cut 30     # tight cut per paper (not the default 200)
  fit.nominal
end

# --- KS0 veto (applied in ROOT analysis, noted here) ---
# The paper vetoes KS0 by doing secondary vertex fit on all π+π- pairs,
# requiring |M(π+π-) - m_KS0| < 0.005 GeV. Events with ≥2 KS0 candidates are rejected.
# In DSL, secondary_vertex_fit on each pair is complex; the veto logic is in ROOT.
alg.note(:KS0_veto, "Secondary vertex fit on all pi+pi- pairs; |M(pipi)-m_KS0|<0.005 GeV; keep events with <2 KS0 candidates")

# --- π0 veto (applied in ROOT analysis) ---
# For events with ≥2 photons, |M(γγ) - m_π0| > 0.01 GeV to veto π0 backgrounds.
alg.note(:pi0_veto, "For events with >=2 photons: |M(gamma gamma) - m_pi0| > 0.01 GeV")

# --- γγ 3(π+π-) background suppression ---
# χ²_4C(J/ψ→γ 6π) < χ²_4C(J/ψ→γγ 6π) comparison applied in ROOT.
alg.note(:gammagamma_suppression, "Require chi2_4C(gamma 6pi) < chi2_4C(gamma gamma 6pi) to suppress gamma gamma 6pi bg")

alg.with_decay_card(<<~DECAYCARD)
  Decay J/psi
  1.0 gamma X1840  PHSP;
  Enddecay
  Decay X1840
  1.0 pi+ pi- pi+ pi- pi+ pi-  PHSP;
  Enddecay
DECAYCARD

alg.apply(sel)
alg.execute_on([data_jpsi])