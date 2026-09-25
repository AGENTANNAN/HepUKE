### Dataset preparation ###
# ψ(3770) at √s = 3.773 GeV — real data and inclusive MC
data_3773  = DatasetManager.real_data.find("712_3773")        # 2.93 fb^-1 ψ(3770) data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")     # matching inclusive MC

# Decay card for the exclusive signal sample:
#   ψ(3770) -> D+ D-, D+ -> ω μ+ ν_μ (ISGW), D- -> K+ π- π-, ω -> π+ π- π0, π0 -> γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-                PHSP;
    Enddecay

    Decay D+
    1.0000 omega mu+ nu_mu      ISGW;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-           PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0          OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# 200,000 exclusive signal events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DptoOmegaMuNu_DmtoKPiPi"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (tag-based, single-tag D-) ###
alg_name = "DpToOmegaMuNuTag"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])          # header file of the algorithm
       .set_constant({"ECMS" => [:double, 3.773]})            # √s = 3.773 GeV
       .with_decay_card(decay_card_signal)

# ---- Tag side: single D- (charm = -1) from the pre-stored DTag candidates ----
tag_alg.tag_side(:Dm) do |t|
  # Six hadronic tag modes: Kππ, K_Sπ, Kπππ0, K_Sππ0, K_Sπππ, KKπ
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                     # pin the tagged side to D-
  # Explicitly requested tag-side windows (store-not-cut is the default)
  t.window :mBC, min: 1.863, max: 1.877          # M_BC in [1.863, 1.877] GeV
  t.window :deltaE, abs: 0.055                   # |ΔE| < 0.055 GeV
end

# ---- Signal side: D+ -> ω μ+ ν_μ (everything the tag did not use) ----
tag_alg.signal_side do |s|
  s.photons 2..2                                  # exactly two good showers (π0 -> γγ)
  s.min_photon_angle 10.0                         # photon–track opening angle > 10°
  s.min_photon_energy 0.025                       # shower energy > 25 MeV
  # μ+ via the ParticleID probability recipe (CL_μ > 0.001, CL_μ > CL_e, CL_μ > CL_K);
  # μ+ , π+ , π-  =>  net signal-side charge +1
  s.charged(mup: 1, pip: 1, pim: 1)
  s.require_charge 1
  s.missing :nu_mu                                # massless missing neutrino
end

# ---- 4C kinematic fit of tag + signal + ν_μ against the measured lab four-momentum ----
tag_alg.fit do |f|
  f.constrain_four_momentum                                                      # 4C constraint
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).between(0.7577, 0.8077)        # |M(π+π-π0) − m_ω| < 0.025
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma, :mup).between(0.0, 1.5)        # M(ω μ+) < 1.5 GeV
  f.chi2_cut 200                                                                 # χ² < 200
end

# ---- BOSS-side steps that cannot be expressed in the tag DSL ----
tag_alg
  .note(:tag_side_deltaE_per_mode, "the π0-containing tag modes (Kπππ0, K_Sππ0, K_Sπππ) use a
    broader |ΔE| window than the 0.055 GeV value; a single tag-side window is declared here
    and the per-mode broadening is applied in the ROOT analysis")
  .note(:background_veto, "signal-side K_S0 veto: reject events with |M(π+π-) − m_K_S0| < 0.015 GeV
    (suppresses D+ -> anti-K0 μ+ ν_μ feed-down)")
  .note(:background_veto, "signal-side extra-photon veto: reject events containing any additional
    photon with energy > 0.15 GeV")
  .note(:background_veto, "signal-side recoil-mass veto: reject events whose recoil mass falls
    inside 0.45–0.55 GeV")
  .note(:pid_correction_method, "muon identification additionally requires EMC energy in
    0.15–0.25 GeV; the lepton PID recipe used by the tag DSL has fixed v1 thresholds and is not
    tunable from the spec")

tag_alg.apply                                   # validate + render (no Selection argument)
root_files = tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])