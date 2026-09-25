# Paper: 2311.05955v3
# Title: First observation of J/ψ → pbar Σ+ KS0 and J/ψ → p Σ- KS0
# Analysis type: ORDINARY
# Dataset: J/ψ (BOSS 708, 3.097 GeV)

# === DATASET PREPARATION ===

data_jpsi = DatasetManager.real_data.find("708_3097")
inc_mc = DatasetManager.load_inclusive_mc.find("708_3097")

# Exclusive MC: J/ψ → pbar Σ+ KS0 + c.c., Σ+ → p π0, π0 → γγ, KS0 → π+π-
exc_mc = DatasetManager.create_exclusive_mc do |c|
  c.decay_card <<~DECAYCARD
    Decay J/psi
    1.0 anti-p- Sigma+ K_S0  PHSP;
    Enddecay
    Decay Sigma+
    1.0 p+ pi0  PHSP;
    Enddecay
    Decay K_S0
    1.0 pi+ pi-  PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma  PHSP;
    Enddecay
  DECAYCARD
  c.related_dataset data_jpsi
end

# === ALGORITHM & SELECTION ===

alg = Algorithm.new("Jpsi_pbarSigma_KS0", version: '00-00-01')
alg.set_header(["Jpsi_pbarSigma_KS0Alg/Jpsi_pbarSigma_KS0.h"])
alg.set_constant({ "ECMS" => [:double, 3.097] })

# Note: The paper reconstructs both charge-conjugate modes:
# (1) J/ψ → pbar Σ+ KS0, Σ+ → p π0, π0 → γγ, KS0 → π+ π-
# (2) J/ψ → p Σ- KS0, Σ- → pbar π0, π0 → γγ, KS0 → π+ π-
# Both share the same Selection chain; yield extraction distinguishes
# pbar Σ+ vs p Σ- in ROOT stage via the proton/antiproton charge.

sel = alg.selection

# --- Track selection ---
sel.select_track do |t|
  t.nTot 4           # 4 charged tracks: p, pbar, π+, π-
  t.nChrp 2          # 2 positive (p from Σ+ or π+ from KS0)
  t.nChrn 2          # 2 negative (pbar or π- from KS0)
end

# --- Photon selection ---
sel.select_photon do |ph|
  ph.nMin 2           # at least 2 photons from π0 → γγ
end

# --- PID: identify proton/antiproton ---
sel.pid(method: :probability) do |pid|
  pid.identify :prp, :prm           # identify one proton, one antiproton
  pid.prob_cut 0                    # loose PID cut per paper
end

# No pion PID applied for KS0 daughters (per paper, Sec IV)

# --- KS0 reconstruction via secondary vertex fit ---
sel.secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) do |ks|
  ks.chi2_max 100
  ks.decay_length_significance 2.0
  # KS0 mass window applied in ROOT; Vz < 20 cm for KS0 daughters
end

# --- π0 reconstruction via Kalman kinematic fit ---
sel.kalman_kinematic_fit([:gamma, :gamma]) do |pi0fit|
  pi0fit.constrain_to :pi0
end

# --- Build Σ+ candidate from proton + π0 ---
sel.build_virtual_particle(:Sigma_plus, from: [:prp, :pi0])

# --- 4C kinematic fit with Σ+ mass constraint ---
# The paper applies a 4C fit plus Σ+ and π0 mass constraints
# In DSL, the kinematic_fit 5C with constrain_mass handles this
sel.kinematic_fit([:prm, :Sigma_plus, :K_S0]) do |fit|
  fit.constrain_four_momentum         # 4C
  fit.constrain_mass(:Sigma_plus)     # Σ+ mass constraint (5C)
  fit.constrain_mass(:pi0)            # π0 mass constraint
  fit.chi2_cut 200                    # loose χ² cut per rule T3
  fit.nominal
end

# --- K_S0 mass selection is in ROOT analysis ---
# The paper extracts KS0 yield from M(π+π-) peak above background

# --- Λ veto: reject events with M(p π-) < 1.126 GeV ---
# (applied in ROOT analysis, not in BOSS selection)

# --- η veto: reject events with M(π+ π- π0) < 0.598 GeV ---
# (applied in ROOT analysis, not in BOSS selection)

alg.with_decay_card(<<~DECAYCARD)
  Decay J/psi
  1.0 anti-p- Sigma+ K_S0  PHSP;
  Enddecay
  Decay Sigma+
  1.0 p+ pi0  PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi-  PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAYCARD

alg.apply(sel)
alg.execute_on([data_jpsi])