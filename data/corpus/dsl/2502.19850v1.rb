# Precision measurement of BF for ψ(2S) → τ+ τ-
# Paper: 2502.19850v1
# Event selection: 2 charged tracks (e + μ), zero photons, missing neutrinos.
# One τ→eνν, the other τ→μνν. Missing mass and cosθ_miss cuts in ROOT.

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  0.1783 mu+ nu_mu anti-nu_tau TAULNUNU;
  0.1783 e+ nu_e anti-nu_tau TAULNUNU;
  Enddecay
  Decay tau-
  0.1783 mu- anti-nu_mu nu_tau TAULNUNU;
  0.1783 e- anti-nu_e nu_tau TAULNUNU;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_tautau"
  config.related_dataset = psip_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

algorithm = Algorithm.new("PsipTauTauAnalysis")
algorithm
  .set_header(["PsipTauTauAnalysisAlg/PsipTauTauAnalysis.h"])
  .set_constant(ECMS: 3.686)

event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==1"
    nChrn       "==1"
    nNet        "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .kinematic_fit([:lp, :lm]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algorithm
  .note(:electron_pid,
    "Electron PID cuts applied in ROOT: E/p(e) ∈ [0.8, 1.2], χ²_dE/dx(e) < 4, |Δtof|(e) < 0.3 ns. " \
    "These are more stringent than the DSL high-momentum lepton defaults.")
  .note(:muon_pid,
    "Muon PID cuts applied in ROOT: E/p(μ) < 0.7, χ²_dE/dx(μ) < 4, |Δtof|(μ) < 0.3 ns. " \
    "Additional MUC depth cut: depth > 81 × (p - 0.65) cm for μ/π separation " \
    "(parameters a=81 cm, b=0.65 GeV/c optimized via FOM = N_sig / sqrt(N_sig + N_bkg)).")
  .note(:track_momentum, "Track momentum < 1.2 GeV/c based on signal MC study.")
  .note(:vertex_fit, "Vertex fit requiring tracks to pass through a common vertex, required to be successful.")
  .note(:missing_mass_cut, "Missing mass M_miss < 3.05 GeV/c² cut applied in ROOT. " \
    "P4_mis = P4_ψ(2S) − P4_eμ, M_miss = sqrt(E_miss² − |p_miss|²).")
  .note(:cos_theta_miss_cut, "|cosθ_miss| < 0.8 cut applied in ROOT. " \
    "cosθ_miss = −(p_e + p_μ)_z / |p_e + p_μ|.")
  .note(:continuum_subtraction,
    "Continuum e+e- → τ+τ- cross section measured from data at √s=3.650 GeV and " \
    "interpolated to ψ(2S) energy. ConExc generator with lineshape used for signal MC.")
  .note(:backgrounds,
    "Main backgrounds: e+e-→e+e-, μ+μ-, γγ (BABAYAGA), two-photon processes (DIAG36, EKHARA, GALUGA). " \
    "QED backgrounds suppressed by lepton PID and missing mass cuts. ψ(2S)→π+π-J/ψ with J/ψ→l+l- " \
    "suppressed by missing mass cut.")
  .note(:signal_mc,
    "Signal MC generated with ConExc using lineshape of Born cross section for e+e- continuum " \
    "in 3.670-3.700 GeV range with beam energy spread effect. τ decays via EVTGEN TAULNUNU model.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)
  .execute_on([psip_data, psip_incMC, exMC_signal])